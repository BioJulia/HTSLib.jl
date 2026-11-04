# The GenomicFeatures extension

With [GenomicFeatures](https://github.com/BioJulia/GenomicFeatures.jl) loaded, alignments and Tabix records have its interval interface, so that code written for any interval can read them:

| Method | Alignment | Tabix record |
|:--|:--|:--|
| `group_name` | `HTSLib.rname(record)` | `HTSLib.sequence_name(record)` |
| `left_position` | `HTSLib.pos(record)` | `record.start + 1` |
| `right_position` | `HTSLib.end_position(record)` | `record.stop` |
| `has_group_name` | whether RNAME isn't `*` | `true` |
| `has_left_position`, `has_right_position` | whether the record has a position | `true` |
| `strand` | `-` with FLAG's 0x10, "SEQ being reverse complemented", `+` without, and `?` if unmapped | `?`, since htslib reads no strand from a line |

Positions are 1-based and inclusive, as GenomicFeatures counts them; [Positions and intervals](@ref) describes the conversion.
A Tabix record always has a group name and both positions, since htslib rejects a line whose sequence the index lacks, or whose columns aren't integers, and reads a start before the first base as the first base.

The interface reads a record as GenomicFeatures reads any interval:

```@example genomicfeatures
using HTSLib
using HTSLib: bgzip, tabix_index, records
import GenomicFeatures
using GenomicFeatures: group_name, left_position, right_position

directory = mktempdir()
path = joinpath(directory, "features.bed")
write(path, "chr1\t1000\t1500\nchr1\t1999\t2100\n")
compressed = bgzip(path)
tabix_index(compressed)

intervals = TabixFile(compressed) do file
    [(group_name(record), left_position(record), right_position(record)) for record in records(file)]
end
```

Views have the interface too, and `group_name` returns the file's `String`, so reading it from a view allocates nothing.

## Conversions to `GenomicInterval`

Records aren't GenomicFeatures' intervals, so its functions that take intervals, such as `is_overlapping`, and its collections, take intervals converted from them.
A conversion is 1-based and inclusive, with the strand above, and keeps the record, unchanged, as the metadata:

| Conversion | Interval | Metadata |
|:--|:--|:--|
| `GenomicInterval(record)` | `GenomicInterval{typeof(record)}`, for an owned record | the record itself |
| `GenomicInterval(view)` | `GenomicInterval{TabixRecord{HTSLib.Owned}}`, or the alignment equivalent | the view's owned copy |
| `GenomicInterval{Nothing}(record)` | `GenomicInterval{Nothing}` | `nothing` |
| `GenomicInterval{T}(record)` | `GenomicInterval{T}` | the record, converted to `T` |

A view's interval holds its copy, since the next iteration overwrites the view's bytes.
The metadata keeps htslib's positions, so a Tabix record's `start` is the interval's `first` less one:

```@example genomicfeatures
using GenomicFeatures: GenomicInterval, GenomicIntervalCollection, is_overlapping, span

interval = TabixFile(compressed) do file
    GenomicInterval(first(records(file)))
end
interval.first, interval.last, interval.metadata.start, interval.metadata.stop
```

A collection takes records, converted to intervals without metadata:

```@example genomicfeatures
collection = TabixFile(compressed) do file
    GenomicIntervalCollection{Nothing}(records(file))
end
length(collection), is_overlapping(interval, GenomicInterval("chr1", 1400, 1600))
```

An alignment without a position has no interval, so converting it throws an `ArgumentError`.
A zero-length BED record, such as `chr1 10 10`, whose start is its stop in htslib, converts to an empty interval, from 11 to 10, whose `span` is 0, and which overlaps nothing.

## Names

HTSLib.jl exports none of its functions, so `using HTSLib, GenomicFeatures` brings no clash: `each_overlap` is GenomicFeatures', and HTSLib.jl's is `HTSLib.each_overlap`, or imported by name.
