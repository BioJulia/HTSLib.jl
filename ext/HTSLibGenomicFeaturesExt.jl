module HTSLibGenomicFeaturesExt

using HTSLib
using GenomicFeatures

# GenomicFeatures positions are 1-based and inclusive.
# Alignments are 1-based already, as SAM's POS is; htslib reports Tabix records 0-based and half-open:
#
# | Record       | GenomicFeatures `group_name` | `left_position` | `right_position` | `strand`                      |
# |:-------------|:-----------------------------|:----------------|:-----------------|:------------------------------|
# | alignment    | `rname`, RNAME               | `pos`, POS      | `end_position`   | from FLAG, or `?` if unmapped |
# | Tabix record | `sequence_name`, `col_seq`   | `start + 1`     | `stop`           | `?`                           |
#
# `group_name` is an alias of the accessor that is true to the file's specification.

GenomicFeatures.group_name(record::AlignmentRecord) = HTSLib.rname(record)
GenomicFeatures.has_group_name(record::AlignmentRecord) = HTSLib.rname(record) !== nothing
GenomicFeatures.left_position(record::AlignmentRecord) = HTSLib.pos(record)
GenomicFeatures.has_left_position(record::AlignmentRecord) = HTSLib.pos(record) > 0
GenomicFeatures.right_position(record::AlignmentRecord) = HTSLib.end_position(record)
GenomicFeatures.has_right_position(record::AlignmentRecord) = HTSLib.pos(record) > 0

# htslib rejects a line whose sequence the index lacks, or whose columns aren't integers, and reads a start before the first base as the first base, so every Tabix record has a group name and both positions.
GenomicFeatures.group_name(record::TabixRecord) = HTSLib.sequence_name(record)
GenomicFeatures.has_group_name(::TabixRecord) = true

GenomicFeatures.left_position(record::TabixRecord) = record.start + 1
GenomicFeatures.has_left_position(::TabixRecord) = true
GenomicFeatures.right_position(record::TabixRecord) = record.stop
GenomicFeatures.has_right_position(::TabixRecord) = true

# The strand of an alignment is that of FLAG, SAM's 0x10, "SEQ being reverse complemented", and unknown for an unmapped record, 0x4.
# htslib doesn't read a strand from a tabix line, whatever its format, so a Tabix record's is unknown.
const FLAG_UNMAPPED = 0x0004
const FLAG_REVERSE = 0x0010

function GenomicFeatures.strand(record::AlignmentRecord)
    flag = HTSLib.core(record).flag
    flag & FLAG_UNMAPPED == 0 || return STRAND_NA
    return flag & FLAG_REVERSE == 0 ? STRAND_POS : STRAND_NEG
end

GenomicFeatures.strand(::TabixRecord) = STRAND_NA

# Conversions to GenomicFeatures' intervals, through its hook, `GenomicInterval{T}(record)`, which converts the record.
# The interval is 1-based and inclusive, from the interval interface above.
# Its metadata is `nothing` for `GenomicInterval{Nothing}`, and otherwise the record, unchanged, with htslib's positions; a view is copied, since the next iteration overwrites its bytes, and an owned record is the metadata itself.
function Base.convert(::Type{GenomicInterval{T}}, record::Union{TabixRecord, AlignmentRecord}) where {T}
    check_position(record)
    return GenomicInterval{T}(group_name(record), left_position(record), right_position(record), strand(record), interval_metadata(T, record))
end

interval_metadata(::Type{Nothing}, record) = nothing
interval_metadata(::Type{T}, record) where {T} = convert(T, copy(record))

check_position(::TabixRecord) = nothing
check_position(record::AlignmentRecord) = has_left_position(record) || throw(ArgumentError("the alignment $(repr(String(HTSLib.qname(record)))) has no position, so it has no interval"))

# `GenomicInterval(record)` keeps an owned record as its metadata, and the owned copy of a view.
GenomicFeatures.GenomicInterval(record::TabixRecord) = GenomicInterval{TabixRecord{HTSLib.Owned}}(record)
GenomicFeatures.GenomicInterval(record::AlignmentRecord) = GenomicInterval{AlignmentRecord{HTSLib.Owned}}(record)

end
