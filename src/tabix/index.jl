# Building Indexes
# ----------------

const PRESETS = (:gff, :bed, :psltbl, :sam, :vcf)

function preset_config(preset::Symbol)
    preset === :gff && return Bindings.tbx_conf_gff()
    preset === :bed && return Bindings.tbx_conf_bed()
    preset === :psltbl && return Bindings.tbx_conf_psltbl()
    preset === :sam && return Bindings.tbx_conf_sam()
    preset === :vcf && return Bindings.tbx_conf_vcf()
    throw(ArgumentError("unknown preset $(repr(preset)); the presets are $(join(map(repr, PRESETS), ", "))"))
end

# The preset that `tabix` picks for a file without `-p`: from its name, ignoring case, or else from its contents if htslib detects VCF.
# Where `tabix` would fall back to GFF, the caller must give a preset.
const PRESET_EXTENSIONS = (".gff.gz" => :gff, ".bed.gz" => :bed, ".sam.gz" => :sam, ".vcf.gz" => :vcf)

# The files that `tabix` detects but indexes with htslib's other indexers, rather than with a tabix preset.
const OTHER_EXTENSIONS = (".bcf" => "BCF", ".bam" => "BAM", ".cram" => "CRAM")

function other_index(path::AbstractString, format::AbstractString)
    indexer = format == "BCF" ? "bcf_index_build" : "bam_index_build"
    return ArgumentError("$(path) is a $(format) file, which htslib indexes with `$(indexer)`, not with a tabix preset; `tabix_index` builds only tabix indexes")
end

function detect_preset(path::AbstractString)
    name = lowercase(path)
    for (extension, preset) in PRESET_EXTENSIONS
        endswith(name, extension) && return preset
    end
    for (extension, format) in OTHER_EXTENSIONS
        endswith(name, extension) && throw(other_index(path, format))
    end
    file = Bindings.hts_open(path, "r")
    if file != C_NULL
        format = unsafe_load(Bindings.hts_get_format(file)).format
        Bindings.hts_close(file)
        format == Bindings.vcf && return :vcf
        format == Bindings.bcf && throw(other_index(path, "BCF"))
        format == Bindings.bam && throw(other_index(path, "BAM"))
        format == Bindings.cram && throw(other_index(path, "CRAM"))
    end
    throw(ArgumentError("cannot tell the preset of $(path) from its name or contents, as `tabix` does; pass `preset`, one of $(join(map(repr, PRESETS), ", "))"))
end

"""
    tabix_index(path::AbstractString; preset = nothing, kwargs...) -> String

Build an index of the BGZF-compressed, sorted file `path` with htslib, as `tabix` does, and return the path of the index.

- `preset`: the columns of the file, as `tabix -p`: `:gff`, `:bed`, `:psltbl`, `:sam` or `:vcf`.
  By default, it is detected as `tabix` detects it: from a name that ends in `.gff.gz`, `.bed.gz`, `.sam.gz` or `.vcf.gz`, ignoring case, or else from the contents of a VCF file.
  Where `tabix` would use the GFF preset for a file that it cannot detect, `tabix_index` throws an `ArgumentError`.
  BAM, CRAM and BCF files, which `tabix` detects but indexes with htslib's other indexers, `bam_index_build` and `bcf_index_build`, throw an `ArgumentError` too, since `tabix_index` builds only tabix indexes.
- `sequence_column`, `begin_column`, `end_column`: 1-based columns that replace those of the preset, as `tabix -s`, `-b` and `-e`.
- `zero_based`: whether positions are 0-based and half-open, instead of 1-based and inclusive; by default, as the preset says.
  `true` is `tabix -0`, and `false` makes a 0-based preset 1-based, as giving `tabix -s`, `-b` and `-e` without `-0` does.
  htslib reads only the begin column differently, since the last base of a 1-based, inclusive interval is the end of the 0-based, half-open one.
- `meta_char`: the character that starts header and comment lines, as `tabix -c`.
- `line_skip`: the number of lines to skip at the start of the file, as `tabix -S`.
- `csi`: build a CSI index, `path * ".csi"`, as `tabix -C`, instead of a TBI index, `path * ".tbi"`.
- `min_shift`: the bins of a CSI index at the finest level are `2^min_shift` bases, as `tabix -m`; 14 by default.
- `index`: the path of the index, instead of the default.
- `threads`: the number of threads that decompress the file, or none.

Whatever the positions of the file, htslib gives each record a 0-based, half-open interval, and [`HTSLib.is_zero_based`](@ref) tells which positions the index was built for.
As htslib does, a start before the first base is read as the first base, and htslib warns "Coordinate <= 0 detected. Did you forget to use the -0 option?", where `-0` is `zero_based = true`.
"""
function tabix_index(
        path::AbstractString;
        preset::Union{Nothing, Symbol} = nothing,
        sequence_column::Union{Nothing, Integer} = nothing,
        begin_column::Union{Nothing, Integer} = nothing,
        end_column::Union{Nothing, Integer} = nothing,
        zero_based::Union{Nothing, Bool} = nothing,
        meta_char::Union{Nothing, Char} = nothing,
        line_skip::Union{Nothing, Integer} = nothing,
        csi::Bool = false,
        min_shift::Integer = csi ? 14 : 0,
        index::Union{Nothing, AbstractString} = nothing,
        threads::Integer = 0,
    )
    csi && min_shift <= 0 && throw(ArgumentError("a CSI index needs a positive min_shift"))
    !csi && min_shift != 0 && throw(ArgumentError("min_shift is for CSI indexes; pass csi = true"))
    config = preset_config(preset === nothing ? detect_preset(path) : preset)
    flags = zero_based === nothing ? config.preset : (zero_based ? config.preset | Bindings.TBX_UCSC : config.preset & ~Bindings.TBX_UCSC)
    config = Bindings.tbx_conf_t(
        flags,
        something(sequence_column, config.sc),
        something(begin_column, config.bc),
        something(end_column, config.ec),
        meta_char === nothing ? config.meta_char : Int32(meta_char),
        something(line_skip, config.line_skip),
    )
    status = Bindings.tbx_index_build3(path, index === nothing ? C_NULL : String(index), min_shift, threads, config)
    status == -2 && error("htslib could not index $(path), since it is not BGZF-compressed")
    status == 0 || error("htslib could not index $(path)")
    return index === nothing ? string(path, csi ? ".csi" : ".tbi") : String(index)
end

# Indexes
# -------

"""
    TabixIndex(path::AbstractString; index = nothing)

The tabix index, TBI or CSI, of the BGZF-compressed file `path`, loaded with htslib.
The index next to the file is found as `tabix` finds it, unless `index` gives its path.
Call `close` when finished; the garbage collector closes it otherwise.
"""
mutable struct TabixIndex
    pointer::Ptr{Bindings.tbx_t}

    function TabixIndex(pointer::Ptr{Bindings.tbx_t})
        index = new(pointer)
        finalizer(close, index)
        return index
    end
end

function TabixIndex(path::AbstractString; index::Union{Nothing, AbstractString} = nothing)
    pointer = Bindings.tbx_index_load3(path, index === nothing ? C_NULL : String(index), 0)
    pointer == C_NULL && error("htslib could not load the index of $(path)" * (index === nothing ? "" : " from $(index)"))
    return TabixIndex(pointer)
end

function Base.close(index::TabixIndex)
    if index.pointer != C_NULL
        Bindings.tbx_destroy(index.pointer)
        index.pointer = C_NULL
    end
    return nothing
end

Base.isopen(index::TabixIndex) = index.pointer != C_NULL

function check_open(index::TabixIndex)
    isopen(index) || throw(ArgumentError("the index is closed"))
    return index.pointer
end

"""
    sequence_names(index::TabixIndex) -> Vector{String}
    sequence_names(file::TabixFile) -> Vector{String}

The sequence names of the index, in the order of the file.
"""
function sequence_names(index::TabixIndex)
    count = Ref{Cint}(0)
    names = Bindings.tbx_seqnames(check_open(index), count)
    names == C_NULL && return String[]
    # The array is allocated by htslib, but the strings belong to the index.
    result = [unsafe_string(unsafe_load(names, k)) for k in 1:count[]]
    Libc.free(names)
    return result
end

# The columns of the indexed file.
config(index::TabixIndex) = unsafe_load(check_open(index)).conf

"""
    HTSLib.is_zero_based(index::TabixIndex) -> Bool
    HTSLib.is_zero_based(file::TabixFile) -> Bool

Whether the positions of the indexed file are 0-based and half-open, as the BED preset and `tabix -0` make them, rather than 1-based and inclusive.
This is htslib's `TBX_UCSC` flag of the index.
Whatever it is, htslib gives each record a 0-based, half-open interval, `start` and `stop`; the line of a record keeps the positions of the file, so this tells how to read its columns.
"""
is_zero_based(index::TabixIndex) = config(index).preset & Bindings.TBX_UCSC != 0
