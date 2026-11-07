# Records and iteration

Every record is one of two kinds, which its type parameter says:
- **an owned record**, `TabixRecord{HTSLib.Owned}` or `AlignmentRecord{HTSLib.Owned}`, whose bytes are in Julia's memory, and which stays valid as long as it's kept, after its file is closed too;
- **a view**, [`TabixRecordView`](@ref) or [`AlignmentRecordView`](@ref), `TabixRecord{HTSLib.Borrowed}` or `AlignmentRecord{HTSLib.Borrowed}`, whose bytes are in a buffer of htslib, which the next record read overwrites, so that it's only valid until the next iteration.

Both are immutable, and a view equals its copy, with the same hash.

## Iterators

`each_record` and `each_overlap` yield owned records by default, each with its bytes copied once; with `copy = false`, they yield views, which allocate nothing:

```@example records
using HTSLib
using HTSLib: bgzip, tabix_index, records, each_record, overlaps, each_overlap

directory = mktempdir()
path = joinpath(directory, "many.bed")
write(path, join(("chr1\t$(10k)\t$(10k + 5)\tfeature$(k)" for k in 1:10_000), '\n'), '\n')
compressed = bgzip(path)
tabix_index(compressed)
file = TabixFile(compressed)

eltype(each_record(file)), eltype(each_record(file; copy = false))
```

A loop over views allocates only to start, whatever the number of records:

```@example records
total_length(records) = sum(record -> record.stop - record.start, records)
total_length(each_record(file; copy = false)) # compile it
@allocated total_length(each_record(file; copy = false))
```

The fields of a view are views too: `HTSLib.line` of a Tabix record is an `ImmutableMemoryView{UInt8}` of htslib's buffer, and `HTSLib.qname` of an alignment a `StringView` of it.
The names of sequences, `HTSLib.sequence_name` and `HTSLib.rname`, are the file's `String`s, shared by its records.
So reading them allocates nothing either.

## Keeping records

A view is overwritten by the next record, so keep one with `copy`, which allocates its bytes once:

```@example records
views = each_record(file; copy = false)
kept = copy(first(views))
next = first(views)
String(kept), String(next)
```

`collect` copies every view, so that `collect(each_record(file; copy = false))` is a vector of owned records, as `records(file)` is.

### Don't keep views

A comprehension, or `push!`, keeps views as they are, so every view kept refers to the same buffer, which holds the last record read.
Copying each view as it comes keeps its own bytes:

```@example records
three_views() = Iterators.take(each_record(file; copy = false), 3)
wrong = [view for view in three_views()]        # views of one buffer
right = [copy(view) for view in three_views()]  # owned records
String.(wrong), String.(right)
```

The three views all read as the third record, while the copies are the three records.
To keep only part of a view, convert it as it comes, as in `[String(view) for view in …]`.

A view of a file that's been closed refers to freed memory, so read it, or copy it, before closing the file.
An owned record keeps what it needs, its line and the file's sequence names, or its data and the file's header.

## Handles

Each iterator opens the file again, and reads through a handle of its own, so iterators of the same file are independent.
Here a query of the first 100 bases, which has 9 records, is interrupted by another query after its first record, and then reads on with its other 8:

```@example records
interrupted = each_overlap(file, "chr1:1-100")
first_record = String(first(interrupted))
other_query = length(overlaps(file, "chr1:50001-50100"))
first_record, other_query, length(collect(interrupted))
```

An iterator closes its handle when it ends, or when it's garbage collected; `close(iterator)` closes it sooner, after a loop that stops early.
Nothing is locked, so an iterator, or a file, is for one task at a time.

```@example records
close(file)
```
