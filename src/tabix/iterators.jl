# Queries
# -------

# Free the buffer that htslib allocated for the lines that an iterator read.
function free_line!(line::Base.RefValue{Bindings.kstring_t})
    Libc.free(line[].s)
    line[] = Bindings.kstring_t()
    return line
end

"""
    each_overlap(file::TabixFile, region::AbstractString; copy = true)
    each_overlap(file::TabixFile, sequence::AbstractString, left::Integer, right::Integer; copy = true)
    each_overlap(file::TabixFile, sequence::AbstractString, range::AbstractUnitRange; copy = true)

Iterate over the records of `file` that overlap a region, as `tabix file region` does, yielding a `TabixRecord{HTSLib.Owned}` of each, which owns its line.
The file must be BGZF-compressed.
The region is a tabix region string, such as `"chr1:1001-2000"` or `"chr1"`, or a sequence name with 1-based, inclusive positions.
As in htslib, a sequence that the index does not have, or a region that htslib cannot parse, has no records.
Each query opens the file again and reads through a handle of its own, so queries of the same file may interleave.

With `copy = false`, it yields a `TabixRecordView` of each, which allocates nothing, and is only valid until the next iteration; `copy` keeps a record.
`collect` gives owned records either way.
"""
function each_overlap(file::TabixFile, region::AbstractString; copy::Bool = true)
    check_bgzf(file)
    return tabix_overlaps(file, Bindings.tbx_itr_querys(file.index.pointer, region), copy)
end

function each_overlap(file::TabixFile, sequence::AbstractString, left::Integer, right::Integer; copy::Bool = true)
    check_bgzf(file)
    tid = Bindings.tbx_name2id(file.index.pointer, sequence)
    tid < 0 && return tabix_overlaps(file, Ptr{Bindings.hts_itr_t}(C_NULL), copy)
    return tabix_overlaps(file, Bindings.tbx_itr_queryi(file.index.pointer, tid, left - 1, right), copy)
end

tabix_overlaps(file, pointer, copy::Bool) = copy ? TabixOverlaps{TabixRecord{Owned}}(file, pointer) : TabixOverlaps{TabixRecordView}(file, pointer)

# Queries seek to the virtual offsets of the index, which only a BGZF-compressed file has.
function check_bgzf(file::TabixFile)
    Bindings.hts_get_bgzfp(check_open(file)) == C_NULL && throw(ArgumentError("queries need a BGZF-compressed file, but $(file.path) is not"))
    return file
end

each_overlap(file::TabixFile, sequence::AbstractString, range::AbstractUnitRange{<:Integer}; copy::Bool = true) = each_overlap(file, sequence, first(range), last(range); copy)

# Iterators read a view of each record, and yield it as a `T`: the view itself, or the record that `copy` gives.
# A query reads the file through a handle of its own, since htslib reads on from the position of the handle within a block of the file, which another query would move.
mutable struct TabixOverlaps{T <: TabixRecord} <: RecordIterator{T}
    const file::TabixFile
    handle::Ptr{Bindings.htsFile}
    pointer::Ptr{Bindings.hts_itr_t}
    const line::Base.RefValue{Bindings.kstring_t}
    const buffer::BorrowedBuffer

    function TabixOverlaps{T}(file::TabixFile, pointer::Ptr{Bindings.hts_itr_t}) where {T}
        handle = Ptr{Bindings.htsFile}(C_NULL)
        if pointer != C_NULL
            handle = try
                open_tabix(file.path, file.threads)
            catch
                Bindings.hts_itr_destroy(pointer)
                rethrow()
            end
        end
        overlaps = new{T}(file, handle, pointer, Ref(Bindings.kstring_t()), BorrowedBuffer())
        finalizer(close, overlaps)
        return overlaps
    end
end

function Base.close(overlaps::TabixOverlaps)
    if overlaps.pointer != C_NULL
        Bindings.hts_itr_destroy(overlaps.pointer)
        overlaps.pointer = C_NULL
    end
    if overlaps.handle != C_NULL
        Bindings.hts_close(overlaps.handle)
        overlaps.handle = C_NULL
    end
    free_line!(overlaps.line)
    return nothing
end

# `iterate` is inlined so that its tuple is built at the call site: returned from a function that is not inlined, a tuple holding a view, which refers to the sequence names, would be allocated on every iteration.
@inline function Base.iterate(overlaps::TabixOverlaps{T}, state = nothing) where {T}
    overlaps.pointer == C_NULL && return nothing
    file = overlaps.file
    # The index belongs to the file, so the file must still be open.
    check_open(file)
    status = Bindings.tbx_itr_next(overlaps.handle, file.index.pointer, overlaps.pointer, overlaps.line)
    if status < 0
        close(overlaps)
        status == -1 && return nothing
        error("htslib could not read a record of $(file.path)")
    end
    current = unsafe_load(overlaps.pointer)
    line = overlaps.line[]
    return convert(T, TabixRecordView(file.sequence_names, current.curr_tid + 1, current.curr_beg, current.curr_end, borrow!(overlaps.buffer, line.s, line.m, line.l))), nothing
end

"""
    overlaps(file::TabixFile, region) -> Vector{TabixRecord{HTSLib.Owned}}
    overlaps(file::TabixFile, sequence, left, right) -> Vector{TabixRecord{HTSLib.Owned}}

The records of `file` that overlap a region, as `collect(each_overlap(file, region...))` gives them.
`query` is an alias of `overlaps`.
"""
overlaps(file::TabixFile, region...) = collect(each_overlap(file, region...))

# Reading Records
# ---------------

"""
    each_record(file::TabixFile; copy = true)

Iterate over every record of `file` from its start, as htslib reads them with `tbx_readrec`, yielding a `TabixRecord{HTSLib.Owned}` of each, which owns its line.
Header lines, the first `line_skip` lines and those that start with the meta character of the index, are skipped, as when the index was built.
The file is opened again, so iterating does not disturb queries.

With `copy = false`, it yields a `TabixRecordView` of each, which allocates nothing, and is only valid until the next iteration; `copy` keeps a record.
`collect` gives owned records either way.
"""
function each_record(file::TabixFile; copy::Bool = true)
    check_open(file)
    return copy ? TabixRecords{TabixRecord{Owned}}(file) : TabixRecords{TabixRecordView}(file)
end

mutable struct TabixRecords{T <: TabixRecord} <: RecordIterator{T}
    const file::TabixFile
    pointer::Ptr{Cvoid}
    line_number::Int
    const line::Base.RefValue{Bindings.kstring_t}
    const tid::Base.RefValue{Cint}
    const start::Base.RefValue{Int64}
    const stop::Base.RefValue{Int64}
    const buffer::BorrowedBuffer

    function TabixRecords{T}(file::TabixFile) where {T}
        pointer = Bindings.bgzf_open(file.path, "r")
        pointer == C_NULL && error("htslib could not open $(file.path)")
        records = new{T}(file, pointer, 0, Ref(Bindings.kstring_t()), Ref{Cint}(0), Ref{Int64}(0), Ref{Int64}(0), BorrowedBuffer())
        finalizer(close, records)
        return records
    end
end

function Base.close(records::TabixRecords)
    if records.pointer != C_NULL
        Bindings.bgzf_close(records.pointer)
        records.pointer = C_NULL
    end
    free_line!(records.line)
    return nothing
end

# Inlined, as `iterate` of `TabixOverlaps` is.
@inline function Base.iterate(records::TabixRecords{T}, state = nothing) where {T}
    records.pointer == C_NULL && return nothing
    file = records.file
    index = check_open(file.index)
    columns = config(file.index)
    while true
        next = Bindings.bgzf_peek(records.pointer)
        if next < 0
            close(records)
            next == -1 && return nothing
            error("htslib could not read $(file.path)")
        end
        if records.line_number < columns.line_skip || next == columns.meta_char
            Bindings.bgzf_getline(records.pointer, Cint('\n'), records.line) < 0 && error("htslib could not read $(file.path)")
            records.line_number += 1
            continue
        end
        status = Bindings.tbx_readrec(records.pointer, index, records.line, records.tid, records.start, records.stop)
        records.line_number += 1
        if status < 0
            close(records)
            status == -1 && return nothing
            error("htslib could not parse line $(records.line_number) of $(file.path), or its sequence isn't in the index")
        end
        line = records.line[]
        return convert(T, TabixRecordView(file.sequence_names, records.tid[] + 1, records.start[], records.stop[], borrow!(records.buffer, line.s, line.m, line.l))), nothing
    end
    return
end

"""
    records(file::TabixFile) -> Vector{TabixRecord{HTSLib.Owned}}

Every record of `file`, as `collect(each_record(file))` gives them.
"""
records(file::TabixFile) = collect(each_record(file))
