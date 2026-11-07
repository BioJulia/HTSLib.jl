# Reading Records
# ---------------

"""
    each_record(file::AlignmentFile; copy = true)

Iterate over every record of `file` from its start, as `samtools view` reads them, yielding an `AlignmentRecord{HTSLib.Owned}` of each, which owns its data.
The file is opened again, so iterating does not disturb queries.

With `copy = false`, it yields an `AlignmentRecordView` of each, which allocates nothing, and is only valid until the next iteration; `copy` keeps a record.
`collect` gives owned records either way.
"""
function each_record(file::AlignmentFile; copy::Bool = true)
    check_open(file)
    return copy ? AlignmentRecords{AlignmentRecord{Owned}}(file) : AlignmentRecords{AlignmentRecordView}(file)
end

# Iterators read a view of each record, and yield it as a `T`: the view itself, or the record that `copy` gives.
mutable struct AlignmentRecords{T <: AlignmentRecord} <: RecordIterator{T}
    const file::AlignmentFile
    pointer::Ptr{Bindings.htsFile}
    header::Ptr{Bindings.sam_hdr_t}
    record::Ptr{Bindings.bam1_t}
    const buffer::BorrowedBuffer

    function AlignmentRecords{T}(file::AlignmentFile) where {T}
        pointer, header = open_alignments(file.path, file.reference)
        alignments = new{T}(file, pointer, header, Bindings.bam_init1(), BorrowedBuffer())
        finalizer(close, alignments)
        return alignments
    end
end

function Base.close(alignments::AlignmentRecords)
    if alignments.record != C_NULL
        Bindings.bam_destroy1(alignments.record)
        alignments.record = C_NULL
    end
    if alignments.header != C_NULL
        Bindings.sam_hdr_destroy(alignments.header)
        alignments.header = C_NULL
    end
    if alignments.pointer != C_NULL
        Bindings.hts_close(alignments.pointer)
        alignments.pointer = C_NULL
    end
    return nothing
end

# `iterate` is inlined so that its tuple is built at the call site: returned from a function that is not inlined, a tuple holding a view, which refers to the header, would be allocated on every iteration.
@inline function Base.iterate(alignments::AlignmentRecords{T}, state = nothing) where {T}
    alignments.pointer == C_NULL && return nothing
    status = Bindings.sam_read1(alignments.pointer, alignments.header, alignments.record)
    if status < 0
        close(alignments)
        status == -1 && return nothing
        error("htslib could not read a record of $(alignments.file.path)")
    end
    return convert(T, borrow_record(alignments.record, alignments.buffer, alignments.file.header)), nothing
end

"""
    records(file::AlignmentFile) -> Vector{AlignmentRecord{HTSLib.Owned}}

Every record of `file`, as `collect(each_record(file))` gives them.
"""
records(file::AlignmentFile) = collect(each_record(file))

# Queries
# -------

"""
    each_overlap(file::AlignmentFile, region::AbstractString; copy = true)
    each_overlap(file::AlignmentFile, sequence::AbstractString, left::Integer, right::Integer; copy = true)
    each_overlap(file::AlignmentFile, sequence::AbstractString, range::AbstractUnitRange; copy = true)

Iterate over the records of `file` that overlap a region, as `samtools view file region` does, yielding an `AlignmentRecord{HTSLib.Owned}` of each, which owns its data.
The region is a samtools region string, such as `"chr1:1001-2000"`, `"chr1"` or `"*"` for the records without a position, or a sequence name with 1-based, inclusive positions.
As in htslib, a sequence that the header does not have has no records.
Queries need the index of the file, and a BAM, CRAM, or SAM file compressed with bgzip.
Each query opens the file again and reads through a handle of its own, so queries of the same file may interleave.

With `copy = false`, it yields an `AlignmentRecordView` of each, which allocates nothing, and is only valid until the next iteration; `copy` keeps a record.
`collect` gives owned records either way.
"""
function each_overlap(file::AlignmentFile, region::AbstractString; copy::Bool = true)
    load_index(file)
    return alignment_overlaps(index -> Bindings.sam_itr_querys(index, file.header.pointer, region), file, copy)
end

function each_overlap(file::AlignmentFile, sequence::AbstractString, left::Integer, right::Integer; copy::Bool = true)
    load_index(file)
    tid = Bindings.sam_hdr_name2tid(file.header.pointer, sequence)
    tid < 0 && return alignment_overlaps(nothing, file, copy)
    return alignment_overlaps(index -> Bindings.sam_itr_queryi(index, tid, left - 1, right), file, copy)
end

each_overlap(file::AlignmentFile, sequence::AbstractString, range::AbstractUnitRange{<:Integer}; copy::Bool = true) = each_overlap(file, sequence, first(range), last(range); copy)

alignment_overlaps(query, file, copy::Bool) = copy ? AlignmentOverlaps{AlignmentRecord{Owned}}(query, file) : AlignmentOverlaps{AlignmentRecordView}(query, file)

# A query reads the file through a handle of its own, with its own header, since htslib reads on from the position of the handle, which another query would move.
# A CRAM index refers to the handle it was loaded for, so a query of a CRAM file loads the index again; the other indexes are the file's.
# `query` makes the iterator from the index, or is `nothing` for a query without records.
mutable struct AlignmentOverlaps{T <: AlignmentRecord} <: RecordIterator{T}
    const file::AlignmentFile
    handle::Ptr{Bindings.htsFile}
    header::Ptr{Bindings.sam_hdr_t}
    # The query's own index, of a CRAM file, or none.
    index::Ptr{Bindings.hts_idx_t}
    pointer::Ptr{Bindings.hts_itr_t}
    record::Ptr{Bindings.bam1_t}
    const buffer::BorrowedBuffer

    function AlignmentOverlaps{T}(query, file::AlignmentFile) where {T}
        overlaps = new{T}(file, C_NULL, C_NULL, C_NULL, C_NULL, Bindings.bam_init1(), BorrowedBuffer())
        finalizer(close, overlaps)
        query === nothing && return overlaps
        try
            overlaps.handle, overlaps.header = open_alignments(file.path, file.reference)
            index = file.index
            if format(file) === :cram
                overlaps.index = read_index(overlaps.handle, file)
                index = overlaps.index
            end
            overlaps.pointer = query(index)
        catch
            close(overlaps)
            rethrow()
        end
        return overlaps
    end
end

function Base.close(overlaps::AlignmentOverlaps)
    if overlaps.record != C_NULL
        Bindings.bam_destroy1(overlaps.record)
        overlaps.record = C_NULL
    end
    if overlaps.pointer != C_NULL
        Bindings.hts_itr_destroy(overlaps.pointer)
        overlaps.pointer = C_NULL
    end
    if overlaps.index != C_NULL
        Bindings.hts_idx_destroy(overlaps.index)
        overlaps.index = C_NULL
    end
    if overlaps.header != C_NULL
        Bindings.sam_hdr_destroy(overlaps.header)
        overlaps.header = C_NULL
    end
    if overlaps.handle != C_NULL
        Bindings.hts_close(overlaps.handle)
        overlaps.handle = C_NULL
    end
    return nothing
end

# Inlined, as `iterate` of `AlignmentRecords` is.
@inline function Base.iterate(overlaps::AlignmentOverlaps{T}, state = nothing) where {T}
    overlaps.pointer == C_NULL && return nothing
    file = overlaps.file
    # The index, but for a CRAM file, belongs to the file, so the file must still be open.
    check_open(file)
    status = Bindings.sam_itr_next(overlaps.handle, overlaps.pointer, overlaps.record)
    if status < 0
        close(overlaps)
        status == -1 && return nothing
        error("htslib could not read a record of $(file.path)")
    end
    return convert(T, borrow_record(overlaps.record, overlaps.buffer, file.header)), nothing
end

"""
    overlaps(file::AlignmentFile, region) -> Vector{AlignmentRecord{HTSLib.Owned}}
    overlaps(file::AlignmentFile, sequence, left, right) -> Vector{AlignmentRecord{HTSLib.Owned}}

The records of `file` that overlap a region, as `collect(each_overlap(file, region...))` gives them.
`query` is an alias of `overlaps`.
"""
overlaps(file::AlignmentFile, region...) = collect(each_overlap(file, region...))
