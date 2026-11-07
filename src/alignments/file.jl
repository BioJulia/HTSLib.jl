# Files
# -----

const FORMAT_NAMES = Dict(Bindings.sam => :sam, Bindings.bam => :bam, Bindings.cram => :cram)

# Open a SAM, BAM or CRAM file and read its header.
function open_alignments(path::AbstractString, reference)
    pointer = Bindings.hts_open(path, "r")
    pointer == C_NULL && error("htslib could not open $(path)")
    try
        format = unsafe_load(Bindings.hts_get_format(pointer)).format
        haskey(FORMAT_NAMES, format) || throw(ArgumentError("$(path) is not a SAM, BAM or CRAM file"))
        reference === nothing || Bindings.hts_set_fai_filename(pointer, reference) == 0 || error("htslib could not load the reference $(reference)")
        header = Bindings.sam_hdr_read(pointer)
        header == C_NULL && error("htslib could not read the header of $(path)")
        return pointer, header
    catch
        Bindings.hts_close(pointer)
        rethrow()
    end
end

"""
    AlignmentFile(path::AbstractString; index = nothing, reference = nothing)
    AlignmentFile(f, path::AbstractString; kwargs...)

An open SAM, BAM or CRAM file, as samtools opens it, with its `header`.

- `index`: the path of the index, BAI, CSI or CRAI, for queries; by default, the index next to the file.
  The index is loaded at the first query.
- `reference`: the path of the FASTA reference of a CRAM file.
  Without it, htslib looks the reference up as samtools does, through the `REF_PATH` and `REF_CACHE` environment variables, which may download it.

As samtools does, htslib adds `MD` and `NM` tags computed from the reference to the records of a CRAM file.

Call `close` when finished, or pass a function `f` to call with the file, which is closed afterwards.
"""
mutable struct AlignmentFile <: AbstractHTSFile
    const path::String
    pointer::Ptr{Bindings.htsFile}
    const header::AlignmentHeader
    index::Ptr{Bindings.hts_idx_t}
    const index_path::Union{Nothing, String}
    const reference::Union{Nothing, String}
end

function AlignmentFile(path::AbstractString; index::Union{Nothing, AbstractString} = nothing, reference::Union{Nothing, AbstractString} = nothing)
    reference = reference === nothing ? nothing : String(reference)
    pointer, header = open_alignments(path, reference)
    index_path = index === nothing ? nothing : String(index)
    file = AlignmentFile(String(path), pointer, AlignmentHeader(header), C_NULL, index_path, reference)
    finalizer(close, file)
    return file
end

function AlignmentFile(f::Function, path::AbstractString; kwargs...)
    file = AlignmentFile(path; kwargs...)
    try
        return f(file)
    finally
        close(file)
    end
end

# The header is left to the garbage collector, since records copied from the file refer to it.
function Base.close(file::AlignmentFile)
    if file.index != C_NULL
        Bindings.hts_idx_destroy(file.index)
        file.index = C_NULL
    end
    if file.pointer != C_NULL
        Bindings.hts_close(file.pointer)
        file.pointer = C_NULL
    end
    return nothing
end

Base.isopen(file::AlignmentFile) = file.pointer != C_NULL

function check_open(file::AlignmentFile)
    isopen(file) || throw(ArgumentError("the file $(file.path) is closed"))
    return file.pointer
end

Base.show(io::IO, file::AlignmentFile) = print(io, "AlignmentFile(", repr(file.path), isopen(file) ? "" : ", closed", ")")

"""
    HTSLib.format(file::AlignmentFile) -> Symbol

The format of the file, `:sam`, `:bam` or `:cram`, as htslib detects it from its content.
"""
format(file::AlignmentFile) = FORMAT_NAMES[unsafe_load(Bindings.hts_get_format(check_open(file))).format]

sequence_names(file::AlignmentFile) = sequence_names(file.header)
sequence_lengths(file::AlignmentFile) = sequence_lengths(file.header)

# Queries need a BGZF-compressed or CRAM file, and its index.
function load_index(file::AlignmentFile)
    pointer = check_open(file)
    flags = unsafe_load(pointer).flags
    is_cram = (flags >> 3) & 0x01 == 0x01
    is_bgzf = (flags >> 4) & 0x01 == 0x01
    (is_bgzf || is_cram) || throw(ArgumentError("queries need a BAM, CRAM or BGZF-compressed SAM file, but $(file.path) is not"))
    file.index == C_NULL && (file.index = read_index(pointer, file))
    return file.index
end

# The index of `file`, loaded for the handle `pointer`, as samtools finds it, or from `file.index_path`.
function read_index(pointer::Ptr{Bindings.htsFile}, file::AlignmentFile)
    index = Bindings.sam_index_load3(pointer, file.path, file.index_path === nothing ? C_NULL : file.index_path, 0)
    index == C_NULL && error("htslib could not load the index of $(file.path)" * (file.index_path === nothing ? "" : " from $(file.index_path)"))
    return index
end
