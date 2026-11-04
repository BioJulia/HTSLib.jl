# Headers
# -------

"""
    AlignmentHeader

The header of a SAM, BAM or CRAM file, as htslib read it.
`String(header)` gives its text, and `sequence_names` and `sequence_lengths` its reference sequences.
It reads the names of its reference sequences once, when it is made, so that every record shares them, and `rname` of a record allocates nothing.
The garbage collector frees it once no file or record refers to it.
"""
mutable struct AlignmentHeader
    pointer::Ptr{Bindings.sam_hdr_t}
    # The names of the reference sequences, by `tid + 1`.
    const sequence_names::Vector{String}

    function AlignmentHeader(pointer::Ptr{Bindings.sam_hdr_t})
        header = new(pointer, read_sequence_names(pointer))
        finalizer(header) do header
            Bindings.sam_hdr_destroy(header.pointer)
            header.pointer = C_NULL
        end
        return header
    end
end

function Base.String(header::AlignmentHeader)
    text = Bindings.sam_hdr_str(header.pointer)
    return text == C_NULL ? "" : unsafe_string(text, Bindings.sam_hdr_length(header.pointer))
end

read_sequence_names(pointer::Ptr{Bindings.sam_hdr_t}) = [unsafe_string(Bindings.sam_hdr_tid2name(pointer, tid)) for tid in 0:(Bindings.sam_hdr_nref(pointer) - 1)]

sequence_names(header::AlignmentHeader) = copy(header.sequence_names)

"""
    sequence_lengths(header::AlignmentHeader) -> Vector{Int}
    sequence_lengths(file::AlignmentFile) -> Vector{Int}

The lengths of the reference sequences, in the order of `sequence_names`.
"""
sequence_lengths(header::AlignmentHeader) = [Int(Bindings.sam_hdr_tid2len(header.pointer, tid)) for tid in 0:(Bindings.sam_hdr_nref(header.pointer) - 1)]

Base.show(io::IO, header::AlignmentHeader) = print(io, "AlignmentHeader(", Bindings.sam_hdr_nref(header.pointer), " reference sequences)")
