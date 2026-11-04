# Generate the functions, structs and constants from the headers of htslib_jll, as `selection.jl` selects them.
# Run with `julia --project=gen gen/generate.jl [--full]`; see README.md.

using Clang
using Clang.Generators
using htslib_jll
using Pkg
using Runic
using TOML

include(joinpath(@__DIR__, "selection.jl"))

const PACKAGE_DIRECTORY = dirname(@__DIR__)
const OUTPUT_DIRECTORY = joinpath(PACKAGE_DIRECTORY, "src", "Bindings")

# Version
# -------

# The bindings are generated from one minor version of htslib, which HTSLib.jl's compat must admit alone.
function check_version()
    version = pkgversion(htslib_jll)
    version = VersionNumber(version.major, version.minor, version.patch)
    compat = TOML.parsefile(joinpath(PACKAGE_DIRECTORY, "Project.toml"))["compat"]["htslib_jll"]
    spec = Pkg.Types.semver_spec(compat)
    version in spec || error("htslib_jll $(version) isn't within HTSLib.jl's compat, \"$(compat)\"; see gen/README.md")
    later = VersionNumber(version.major, version.minor + 1)
    earlier = VersionNumber(version.major, version.minor - 1, typemax(Int32))
    (later in spec || (version.minor > 0 && earlier in spec)) &&
        error("HTSLib.jl's compat, \"$(compat)\", admits other minor versions than htslib_jll $(version), which the bindings are generated from; see gen/README.md")
    return
end

# Parsing
# -------

function parse_headers()
    include_directory = joinpath(htslib_jll.artifact_dir, "include")
    # bgzf.h defines `struct BGZF`, which hts.h only declares, so it comes first, or Clang.jl's sanity check fails.
    headers = [joinpath(include_directory, "htslib", name) for name in ("bgzf.h", "hts.h", "sam.h", "tbx.h")]
    options = load_options(joinpath(@__DIR__, "generator.toml"))
    options["general"]["output_file_path"] = joinpath(@__DIR__, "Generated.jl")
    arguments = get_default_args()
    push!(arguments, "-I$(include_directory)")
    return Generators.create_context(headers, arguments, options)
end

# The nodes of each identifier, of which there may be several, such as a struct and its typedef.
function index_nodes(context)
    nodes = Dict{Symbol, Vector{Generators.ExprNode}}()
    for node in context.dag.nodes
        push!(get!(nodes, node.id, Generators.ExprNode[]), node)
    end
    return nodes
end

function find_node(predicate, nodes, id)
    for node in get(nodes, id, Generators.ExprNode[])
        predicate(node) && return node
    end
    return error("htslib has no $(id) of the kind that selection.jl expects")
end

is_function(node) = node.type isa Generators.AbstractFunctionNodeType && !isempty(node.exprs)
is_struct(node) = node.type isa Generators.AbstractStructNodeType && !isempty(node.exprs)
is_macro(node) = node.type isa Generators.MacroDefault
is_enum(node) = node.type isa Generators.AbstractEnumNodeType

# The values of all enum members, by name.
function enum_values(context)
    values = Dict{Symbol, Int}()
    for node in context.dag.nodes
        is_enum(node) || continue
        for child in Clang.children(node.cursor)
            child isa CLEnumConstantDecl && (values[Symbol(Clang.spelling(child))] = Clang.value(child))
        end
    end
    return values
end

enum_names(context) = Set(node.id for node in context.dag.nodes if is_enum(node))

# Functions
# ---------

# Identifiers in front of declarations that expand to attributes.
const ANNOTATIONS = Set(["HTSLIB_EXPORT", "HTS_RESULT_USED", "HTS_DEPRECATED"])

# The C declaration of a function, as written in its header, on one line.
function declaration(cursor)
    tokens = [token.text for token in Clang.tokenize(cursor) if token.text ∉ ANNOTATIONS]
    text = join(tokens, " ")
    text = replace(text, r"\s*\*\s*" => " *", r"\s*\(\s*" => "(", r"\s*\)" => ")", r"\s*,\s*" => ", ", r"\s*\[\s*" => "[", r"\s*\]" => "]")
    text = replace(text, "* *" => "**")
    return text * ";"
end

# The `@ccall` expression in a function that Clang.jl generated.
function find_ccall(expression)
    expression isa Expr || return nothing
    expression.head === :macrocall && expression.args[1] === Symbol("@ccall") && return expression
    for argument in expression.args
        found = find_ccall(argument)
        found === nothing || return found
    end
    return nothing
end

function rewrite_type(type, rules, enums)
    type isa Symbol && type in enums && return :Cint
    return get(rules, type, type)
end

function emit_function(io, nodes, enums, name, overrides)
    node = find_node(is_function, nodes, name)
    call = find_ccall(node.exprs[1])
    signature = call.args[end]
    returns = signature.args[2]
    arguments = signature.args[1].args[2:end]

    argument_overrides = get(overrides, :arguments, (;))
    for key in keys(argument_overrides)
        any(argument -> argument.args[1] === key, arguments) || error("$(name) has no argument $(key)")
    end

    names = Symbol[]
    typed = String[]
    for argument in arguments
        c_name, type = argument.args
        type = get(argument_overrides, c_name, rewrite_type(type, ARGUMENT_TYPES, enums))
        julia_name = get(ARGUMENT_NAMES, c_name, c_name)
        push!(names, julia_name)
        push!(typed, "$(julia_name)::$(type)")
    end
    returns = get(overrides, :returns, rewrite_type(returns, RETURN_TYPES, enums))

    call_text = "@ccall libhts.$(name)($(join(typed, ", ")))::$(returns)"
    haskey(overrides, :wrap) && (call_text = "$(overrides.wrap)($(call_text))")
    comment = declaration(node.cursor)
    haskey(overrides, :note) && (comment *= "  // $(overrides.note)")

    println(io, "# ", comment)
    definition = "$(name)($(join(names, ", "))) = $(call_text)"
    if length(definition) <= 140
        println(io, definition)
    else
        println(io, "function $(name)($(join(names, ", ")))")
        println(io, "    return $(call_text)")
        println(io, "end")
    end
    println(io)
    return
end

entry(name::Symbol) = (name, (;))
entry(pair::Pair) = (pair.first, pair.second)

# Structs and constants
# ---------------------

function emit_struct(io, nodes, name, overrides)
    node = find_node(is_struct, nodes, name)
    expression = Base.remove_linenums!(deepcopy(node.exprs[1]))
    expression.head === :struct || error("$(name) isn't a plain struct; describe it in LAYOUTS instead")
    field_overrides = get(overrides, :fields, (;))
    fields = expression.args[3].args
    for key in keys(field_overrides)
        any(field -> field.args[1] === key, fields) || error("$(name) has no field $(key)")
    end
    println(io, "\"\"\"")
    println(io, "    ", name)
    println(io)
    println(io, overrides.doc)
    println(io, "\"\"\"")
    println(io, "struct ", name)
    for field in fields
        field_name, type = field.args
        println(io, "    ", field_name, "::", get(field_overrides, field_name, type))
    end
    println(io, "end")
    println(io)
    return
end

function emit_constant(io, nodes, enums, name, type)
    value = if haskey(enums, name)
        enums[name]
    else
        node = find_node(is_macro, nodes, name)
        node.exprs[1].args[1].args[2]
    end
    getfield(Base, type)(value)  # errors if the value doesn't fit the type
    println(io, "const ", name, " = ", type, "(", repr(value), ")")
    return
end

# Layouts
# -------

# The offset in bytes of a field, given by a path such as `"version.minor"` through nested structs, or -1 if there's no such field.
function c_offset(c_type, path)
    offset = 0
    for name in split(path, '.')
        bits = Clang.getOffsetOf(c_type, name)
        bits < 0 && return -1
        offset += bits ÷ 8
        declaration = Clang.getTypeDeclaration(c_type)
        index = findfirst(child -> child isa CLFieldDecl && Clang.spelling(child) == name, Clang.children(declaration))
        c_type = Clang.getCanonicalType(Clang.getCursorType(Clang.children(declaration)[index]))
    end
    return offset
end

# Compare the hand-written structs of Bindings with the offsets and sizes that libclang computes.
function check_layouts(nodes)
    sandbox = Module(:LayoutCheck)
    Base.include(sandbox, joinpath(PACKAGE_DIRECTORY, "src", "Bindings.jl"))
    Bindings = Base.invokelatest(getfield, sandbox, :Bindings)
    differences = String[]
    for (name, layout) in LAYOUTS
        node = find_node(node -> node.type isa Generators.AbstractStructNodeType && !(node.type isa Generators.StructForwardDecl), nodes, name)
        c_type = Clang.getCanonicalType(Clang.getCursorType(node.cursor))
        julia_type = Base.invokelatest(getfield, Bindings, name)
        for (index, field) in enumerate(fieldnames(julia_type))
            c_field = string(get(layout.fields, field, field))
            offset = c_offset(c_type, c_field)
            if offset < 0
                push!(differences, "$(name).$(field): htslib has no field $(c_field)")
            elseif offset != fieldoffset(julia_type, index)
                push!(differences, "$(name).$(field): offset $(fieldoffset(julia_type, index)) in Bindings, $(offset) in htslib")
            end
        end
        if !layout.prefix && sizeof(julia_type) != Clang.getSizeOf(c_type)
            push!(differences, "$(name): size $(sizeof(julia_type)) in Bindings, $(Clang.getSizeOf(c_type)) in htslib")
        end
    end
    isempty(differences) || error("The layouts of Bindings differ from htslib's:\n" * join(differences, "\n"))
    @info "Layouts of $(join(first.(LAYOUTS), ", ")) agree with htslib"
    return
end

# Output
# ------

function header()
    return """
    # This file was generated by gen/generate.jl from htslib_jll $(pkgversion(htslib_jll)), of what gen/selection.jl selects, and is included by src/Bindings.jl.
    # Don't edit it; edit gen/selection.jl and run the generator instead.

    """
end

function section(io, title)
    println(io, "# ", title)
    println(io, "# ", "-"^length(title))
    println(io)
    return
end

function write_formatted(path, text)
    write(path, Runic.format_string(rstrip(text) * "\n"))
    @info "Wrote $(relpath(path, PACKAGE_DIRECTORY))"
    return
end

function main(arguments)
    check_version()
    context = parse_headers()
    build!(context, BUILDSTAGE_NO_PRINTING)
    nodes = index_nodes(context)
    enums = enum_names(context)
    values = enum_values(context)

    types = IOBuffer()
    print(types, header())
    for (name, overrides) in STRUCTS
        emit_struct(types, nodes, name, overrides)
    end
    for (comment, constants) in CONSTANTS
        println(types, "# ", comment)
        for (name, type) in constants
            emit_constant(types, nodes, values, name, type)
        end
        println(types)
    end

    functions = IOBuffer()
    print(functions, header())
    println(functions, "# Optional strings are `Ptr{UInt8}`, so that `C_NULL` can be passed; a `String` converts to a pointer to its NUL-terminated bytes.")
    println(functions)
    for (title, entries) in FUNCTIONS
        section(functions, title)
        for item in entries
            emit_function(functions, nodes, enums, entry(item)...)
        end
    end

    mkpath(OUTPUT_DIRECTORY)
    write_formatted(joinpath(OUTPUT_DIRECTORY, "types.jl"), String(take!(types)))
    write_formatted(joinpath(OUTPUT_DIRECTORY, "functions.jl"), String(take!(functions)))
    check_layouts(nodes)

    "--full" in arguments && build!(context, BUILDSTAGE_PRINTING_ONLY)
    return
end

main(ARGS)
