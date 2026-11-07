# Files
# -----

"""
    TabixFile(path::AbstractString; index = nothing, threads = 0)
    TabixFile(f, path::AbstractString; kwargs...)

An open file, BGZF-compressed or plain, with its tabix index, as `tabix` opens it.
The index is found next to the file, unless `index` gives its path.
A plain file needs `index`, such as that of its compressed copy, and can be read from start to end with `each_record`, but not queried.
`threads` is the number of threads that decompress the file, or none, for the file and for each query.

Each query, and each read from start to end, opens the file again and reads through a handle of its own, so that they may interleave.

Call `close` when finished, or pass a function `f` to call with the file, which is closed afterwards.
"""
mutable struct TabixFile <: AbstractHTSFile
    const path::String
    pointer::Ptr{Cvoid}
    const index::TabixIndex
    const sequence_names::Vector{String}
    const threads::Int
end

# Open `path` with htslib, with `threads` threads that decompress it, or none.
function open_tabix(path::AbstractString, threads::Integer)
    pointer = Bindings.hts_open(path, "r")
    pointer == C_NULL && error("htslib could not open $(path)")
    if threads > 0 && Bindings.hts_set_threads(pointer, threads) != 0
        Bindings.hts_close(pointer)
        error("htslib could not start $(threads) threads to read $(path)")
    end
    return pointer
end

function TabixFile(path::AbstractString; index::Union{Nothing, AbstractString} = nothing, threads::Integer = 0)
    loaded = TabixIndex(path; index)
    pointer = try
        open_tabix(path, threads)
    catch
        close(loaded)
        rethrow()
    end
    file = TabixFile(String(path), pointer, loaded, sequence_names(loaded), Int(threads))
    finalizer(close, file)
    return file
end

function TabixFile(f::Function, path::AbstractString; kwargs...)
    file = TabixFile(path; kwargs...)
    try
        return f(file)
    finally
        close(file)
    end
end

function Base.close(file::TabixFile)
    if file.pointer != C_NULL
        Bindings.hts_close(file.pointer)
        file.pointer = C_NULL
    end
    close(file.index)
    return nothing
end

Base.isopen(file::TabixFile) = file.pointer != C_NULL

function check_open(file::TabixFile)
    isopen(file) || throw(ArgumentError("the file $(file.path) is closed"))
    return file.pointer
end

sequence_names(file::TabixFile) = copy(file.sequence_names)
is_zero_based(file::TabixFile) = is_zero_based(file.index)

Base.show(io::IO, file::TabixFile) = print(io, "TabixFile(", repr(file.path), isopen(file) ? "" : ", closed", ")")
