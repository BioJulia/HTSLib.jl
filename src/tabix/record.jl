# Records
# -------

"""
    TabixRecord{O}

A record that htslib read from a tabix-indexed file: the 1-based index of its sequence in `sequence_names(file)`, its interval, and its line.
`O` is its ownership:
- a `TabixRecordView`, which is `TabixRecord{HTSLib.Borrowed}`, borrows a buffer of htslib that the next record overwrites, so it is only valid until the next iteration;
- a `TabixRecord{HTSLib.Owned}` owns its line, in Julia's memory, so it stays valid after its file is closed.

`copy(view)` gives an owned record, which allocates only its line, and since an owned record is immutable, `copy(record)` is `record` itself.
A view equals its copy.
`TabixRecord(sequence, start, stop, line)` makes an owned record from strings.

The fields `start` and `stop` are the interval, 0-based and half-open, as htslib reports it: `start` is the position of the first base, and `stop` the position after the last.
`sequence_name(record)` gives its sequence name, from the column that the tabix specification calls `col_seq`, the "Column for the sequence name".
`position(record)` gives the 1-based position of the first base, and `interval_0based(record)` the interval.
`line(record)` gives the line as an `ImmutableMemoryView{UInt8}`, without its line ending, and `String(record)` copies it.
The line keeps the positions of the file, which `HTSLib.is_zero_based(file)` tells.
With GenomicFeatures loaded, `group_name`, `left_position` and `right_position` give the sequence name and the 1-based, inclusive interval.
"""
struct TabixRecord{O <: Ownership} <: AbstractHTSRecord
    # The file's sequence names, shared by its records rather than copied for each.
    sequence_names::Vector{String}
    sequence_index::Int
    start::Int
    stop::Int
    line::ImmutableMemoryView{UInt8}
end

"""
    TabixRecordView

A [`TabixRecord`](@ref) that borrows a buffer of htslib, `TabixRecord{HTSLib.Borrowed}`, which is only valid until the next iteration.
"""
const TabixRecordView = TabixRecord{Borrowed}

TabixRecord(sequence::AbstractString, start::Integer, stop::Integer, line::AbstractString) = TabixRecord{Owned}([String(sequence)], 1, start, stop, ImmutableMemoryView(String(line)))

Base.copy(record::TabixRecordView) = TabixRecord{Owned}(record.sequence_names, record.sequence_index, record.start, record.stop, copy(record.line))
Base.copy(record::TabixRecord{Owned}) = record

# Iterators yield a view as the type of their elements, which for an owned record copies it.
Base.convert(::Type{TabixRecord{Owned}}, record::TabixRecordView) = copy(record)
owned(::Type{<:TabixRecord}) = TabixRecord{Owned}

"""
    HTSLib.sequence_name(record::TabixRecord) -> String

The sequence name of `record`, from the column that the tabix specification calls `col_seq`, the "Column for the sequence name".
It is the file's `String`, shared by its records, so reading it allocates nothing.
"""
sequence_name(record::TabixRecord) = record.sequence_names[record.sequence_index]

"""
    HTSLib.position(record::TabixRecord) -> Int

The 1-based position of the first base of `record`, `record.start + 1`, whatever the positions of the file.
"""
position(record::TabixRecord) = record.start + 1

"""
    HTSLib.interval_0based(record::TabixRecord) -> Tuple{Int, Int}

The interval of `record` as htslib gives it, 0-based and half-open, `(record.start, record.stop)`, whatever the positions of the file.
"""
interval_0based(record::TabixRecord) = (record.start, record.stop)
"""
    HTSLib.line(record::TabixRecord) -> ImmutableMemoryView{UInt8}

The line of `record`, without its line ending, as bytes that can't be written through the view.
A view's line is in a buffer of htslib, so it is only valid until the next iteration; `String(record)` copies it.
`split_each(line(record), UInt8('\\t'))`, of MemoryViews, gives its columns without copying them, and `StringView(line(record))`, of StringViews, gives it as a string without copying it.
"""
line(record::TabixRecord) = record.line
Base.String(record::TabixRecord) = String(record.line)

Base.:(==)(a::TabixRecord, b::TabixRecord) = a.start == b.start && a.stop == b.stop && sequence_name(a) == sequence_name(b) && a.line == b.line
Base.hash(record::TabixRecord, h::UInt) = hash(record.line, hash(record.stop, hash(record.start, hash(sequence_name(record), hash(TabixRecord, h)))))

record_name(::TabixRecordView) = "TabixRecordView"
record_name(::TabixRecord{Owned}) = "TabixRecord"
# The line is shown through a `StringView`, without copying it.
function Base.show(io::IO, record::TabixRecord)
    print(io, record_name(record), "(")
    show(io, sequence_name(record))
    print(io, ", ", record.start, ", ", record.stop, ", ")
    show(io, StringView(record.line))
    return print(io, ")")
end
