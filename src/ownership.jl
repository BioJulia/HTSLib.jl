# Ownership
# ---------
#
# A record either owns its bytes, in Julia's memory, or borrows a buffer of htslib, which the next record read overwrites.
# The two differ only in who owns the bytes, so they are one type, with the owner as a type parameter, as `MemoryView` has its mutability.

abstract type Ownership end

"""
    HTSLib.Owned

The type parameter of a record that owns its bytes, in Julia's memory, so that it stays valid after the next iteration and after its file is closed: `TabixRecord{Owned}`.
`copy` gives one from a view.
"""
struct Owned <: Ownership end

"""
    HTSLib.Borrowed

The type parameter of a record that borrows a buffer of htslib, which the next record read overwrites, so that it is only valid until the next iteration: `TabixRecordView`, which is `TabixRecord{Borrowed}`.
"""
struct Borrowed <: Ownership end

# A buffer of htslib, as `Memory` that wraps it without owning it.
# It is wrapped again only when htslib reallocates the buffer, so that the views of each record allocate nothing.
mutable struct BorrowedBuffer
    pointer::Ptr{UInt8}
    memory::Memory{UInt8}
end

BorrowedBuffer() = BorrowedBuffer(C_NULL, Memory{UInt8}())

# A view of the first `n` bytes of the buffer at `pointer`, which holds `capacity` bytes.
function borrow!(buffer::BorrowedBuffer, pointer::Ptr{UInt8}, capacity::Integer, n::Integer)
    capacity = max(Int(capacity), Int(n))
    if pointer != buffer.pointer || capacity != length(buffer.memory)
        buffer.memory = pointer == C_NULL ? Memory{UInt8}() : unsafe_wrap(Memory{UInt8}, pointer, capacity)
        buffer.pointer = pointer
    end
    return ImmutableMemoryView(buffer.memory)[1:Int(n)]
end

# Iterators
# ---------

# The type of a record that owns its bytes, of the same kind as the record type `T`.
function owned end

# An iterator over the records of a file, which yields a `T`: an owned record, or a view.
# Collecting views would keep views of one buffer, so `collect` copies them, as owned records.
abstract type RecordIterator{T <: AbstractHTSRecord} end

Base.IteratorSize(::Type{<:RecordIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:RecordIterator{T}}) where {T} = T
Base.collect(iterator::RecordIterator{T}) where {T} = collect(owned(T), iterator)
