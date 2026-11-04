# Records
# -------

"""
    AlignmentRecord{O}

An alignment that htslib read: the fixed-length fields of its `bam1_core_t`, its variable-length data, and the header of its file.
`O` is its ownership:
- an `AlignmentRecordView`, which is `AlignmentRecord{HTSLib.Borrowed}`, borrows a buffer of htslib that the next record overwrites, so it is only valid until the next iteration;
- an `AlignmentRecord{HTSLib.Owned}` owns its data, in Julia's memory, so it stays valid after its file is closed.

`copy(view)` gives an owned record, which allocates only its data, and since an owned record is immutable, `copy(record)` is `record` itself.
Two records are equal if they have the same fields and come from the same header, so a view equals its copy.

`string(record)` is the SAM line of the record, as `samtools view` prints it, with every field.
`data(record)` gives the variable-length fields, QNAME, CIGAR, SEQ, QUAL and the tags, as htslib lays them out, as an `ImmutableMemoryView{UInt8}`.
The accessors take the names of the SAM specification's fields, with aliases where the specification names a field otherwise:

| Accessor | Value | Aliases |
|:--|:--|:--|
| `qname` | QNAME, "Query template NAME", as a `StringView` of the record's data, valid as long as the record is | `template_name` |
| `rname` | RNAME, "Reference sequence NAME", as the header's `String`, or `nothing` without one | `reference_name` |
| `pos` | POS, "1-based leftmost mapping POSition", or 0 without one | `position` |
| `end_position` | the last reference base that the CIGAR covers, 1-based, as htslib's `bam_endpos` computes it: an unmapped record, or one without a CIGAR, covers one base, and a record without a position ends at 0; derived, not a field | — |
| `reference_index` | the 1-based index of RNAME into `sequence_names`, or 0 without one; derived, not a field | — |

A record overlaps a region from `left` to `right` if `pos(record) <= right` and `end_position(record) >= left`.
With GenomicFeatures loaded, `group_name`, `left_position` and `right_position` give `rname`, `pos` and `end_position`.
"""
struct AlignmentRecord{O <: Ownership} <: AbstractHTSRecord
    core::Bindings.bam1_core_t
    # The variable-length fields, QNAME, CIGAR, SEQ, QUAL and the tags, as htslib lays them out.
    data::ImmutableMemoryView{UInt8}
    header::AlignmentHeader
end

"""
    AlignmentRecordView

An [`AlignmentRecord`](@ref) that borrows a buffer of htslib, `AlignmentRecord{HTSLib.Borrowed}`, which is only valid until the next iteration.
"""
const AlignmentRecordView = AlignmentRecord{Borrowed}

# A view of the record that htslib read into `bam`, whose data is in `buffer`.
function borrow_record(bam::Ptr{Bindings.bam1_t}, buffer::BorrowedBuffer, header::AlignmentHeader)
    record = unsafe_load(bam)
    return AlignmentRecordView(record.core, borrow!(buffer, record.data, record.m_data, record.l_data), header)
end

Base.copy(record::AlignmentRecordView) = AlignmentRecord{Owned}(record.core, copy(record.data), record.header)
Base.copy(record::AlignmentRecord{Owned}) = record

# Iterators yield a view as the type of their elements, which for an owned record copies it.
Base.convert(::Type{AlignmentRecord{Owned}}, record::AlignmentRecordView) = copy(record)
owned(::Type{<:AlignmentRecord}) = AlignmentRecord{Owned}

Base.:(==)(a::AlignmentRecord, b::AlignmentRecord) = a.header === b.header && a.core === b.core && a.data == b.data
Base.hash(record::AlignmentRecord, h::UInt) = hash(record.data, hash(record.core, hash(objectid(record.header), hash(AlignmentRecord, h))))

# htslib's `BAM_USER_OWNS_DATA`, of its `sam.h`: the caller owns the data of a `bam1_t`, which htslib then neither frees nor reallocates.
const BAM_USER_OWNS_DATA = UInt32(2)

# Call `f` with a pointer to a `bam1_t` of `record`, which points at its data, and lives only during the call.
function with_bam(f, record::AlignmentRecord)
    data = record.data
    GC.@preserve data begin
        bam = Ref(Bindings.bam1_t(record.core, 0, pointer(data), length(data), length(data), BAM_USER_OWNS_DATA))
        return GC.@preserve bam f(Base.unsafe_convert(Ptr{Bindings.bam1_t}, bam))
    end
end

core(record::AlignmentRecord) = record.core
"""
    HTSLib.data(record::AlignmentRecord) -> ImmutableMemoryView{UInt8}

The variable-length fields of `record`, QNAME, CIGAR, SEQ, QUAL and the tags, as htslib lays them out in its `bam1_t`, as bytes that can't be written through the view.
A view's data is in a buffer of htslib, so it is only valid until the next iteration; `copy(record)` keeps it.
"""
data(record::AlignmentRecord) = record.data

# The fields of the SAM specification.
# QNAME comes first in the data, with a NUL and up to three more, `l_extranul`, which align what follows.

"""
    HTSLib.qname(record::AlignmentRecord) -> StringView

QNAME, the SAM specification's "Query template NAME", as a `StringView` of the record's data, so that reading it allocates nothing.
It is valid as long as the record is: until the next iteration for a view; `String(qname(record))` keeps it.
`template_name` is an alias.
"""
qname(record::AlignmentRecord) = StringView(record.data[1:(Int(record.core.l_qname) - Int(record.core.l_extranul) - 1)])

"""
    HTSLib.pos(record::AlignmentRecord) -> Int

POS, the SAM specification's "1-based leftmost mapping POSition" of the first base that the CIGAR covers, or 0 for a record without a position.
`position` is an alias.
"""
pos(record::AlignmentRecord) = Int(core(record).pos) + 1

"""
    HTSLib.rname(record::AlignmentRecord) -> Union{String, Nothing}

RNAME, the SAM specification's "Reference sequence NAME", or `nothing` for a record without one, whose RNAME is `*`.
It is the header's `String`, shared by the records of the file, so reading it allocates nothing.
`reference_name` is an alias.
"""
function rname(record::AlignmentRecord)
    tid = core(record).tid
    return tid < 0 ? nothing : record.header.sequence_names[tid + 1]
end

# Aliases from the SAM specification's descriptions of the fields: "Query template NAME", "Reference sequence NAME" and "mapping POSition".

"""
    HTSLib.template_name(record::AlignmentRecord) -> StringView

An alias of [`HTSLib.qname`](@ref), after the SAM specification's "Query template NAME".
"""
template_name(record::AlignmentRecord) = qname(record)

"""
    HTSLib.reference_name(record::AlignmentRecord) -> Union{String, Nothing}

An alias of [`HTSLib.rname`](@ref), after the SAM specification's "Reference sequence NAME".
"""
reference_name(record::AlignmentRecord) = rname(record)

"""
    HTSLib.position(record::AlignmentRecord) -> Int

An alias of [`HTSLib.pos`](@ref), after the SAM specification's "mapping POSition".
"""
position(record::AlignmentRecord) = pos(record)

# Derived values, which are not fields.

"""
    HTSLib.reference_index(record::AlignmentRecord) -> Int

The 1-based index of RNAME into `sequence_names` of the file, or 0 for a record without one.
"""
reference_index(record::AlignmentRecord) = Int(core(record).tid) + 1

"""
    HTSLib.end_position(record::AlignmentRecord) -> Int

The 1-based position of the last reference base that the CIGAR covers, as htslib's `bam_endpos` computes it.
As htslib does, an unmapped record, or one without a CIGAR, covers one base, `pos(record)`, and a record without a position ends at 0.
"""
end_position(record::AlignmentRecord) = Int(with_bam(Bindings.bam_endpos, record))

"""
    HTSLib.interval_0based(record::AlignmentRecord) -> Tuple{Int, Int}

The interval of `record` as htslib gives it, 0-based and half-open: `(pos(record) - 1, end_position(record))`, which is `(-1, 0)` for a record without a position, as htslib indexes it.
"""
interval_0based(record::AlignmentRecord) = (pos(record) - 1, end_position(record))

function Base.string(record::AlignmentRecord)
    line = Ref(Bindings.kstring_t())
    try
        with_bam(bam -> Bindings.sam_format1(record.header.pointer, bam, line), record) < 0 && error("htslib could not format a record")
        return unsafe_string(line[].s, line[].l)
    finally
        Libc.free(line[].s)
    end
end

record_name(::AlignmentRecordView) = "AlignmentRecordView"
record_name(::AlignmentRecord{Owned}) = "AlignmentRecord"
Base.show(io::IO, record::AlignmentRecord) = print(io, record_name(record), "(", repr(string(record)), ")")
