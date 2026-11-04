# Positions and intervals

Genomic files count positions in two ways:
- **1-based and inclusive**, as GFF, VCF and SAM do, and as people usually write a region: the bases 1001 to 2000 of `chr1` are `chr1:1001-2000`;
- **0-based and half-open**, as BED does, and as htslib does internally: the same bases are the interval from 1000 to 2000, which includes the base at 1000, the 1001st, and excludes the one at 2000.

The two describe the same bases with the same end, and starts that differ by one: an interval from `start` to `stop`, 0-based and half-open, is the bases `start + 1` to `stop`, 1-based and inclusive.
HTSLib.jl gives positions as htslib does, so which one a value uses depends on where it comes from.

## Where each is used

| Value | Positions |
|:--|:--|
| `record.start` and `record.stop` of a [`TabixRecord`](@ref) | 0-based, half-open, as htslib gives them, whatever the file's |
| [`HTSLib.position`](@ref) of a Tabix record | 1-based, `start + 1` |
| [`HTSLib.interval_0based`](@ref) of a Tabix record | 0-based, half-open, `(start, stop)` |
| the columns of `HTSLib.line(record)` | the file's own, which [`HTSLib.is_zero_based`](@ref) tells |
| region strings, and the positions of a query | 1-based, inclusive |
| [`HTSLib.pos`](@ref) of an alignment | 1-based, SAM's POS |
| [`HTSLib.end_position`](@ref) of an alignment | 1-based, the last base |
| [`HTSLib.interval_0based`](@ref) of an alignment | 0-based, half-open, `(pos - 1, end_position)` |
| `left_position` and `right_position`, of the GenomicFeatures extension | 1-based, inclusive |

## How htslib reads a tabix line

For each line, htslib reads the sequence name, the start and the end from the columns of the index:
- **the start:** the begin column, less one for a 1-based file, so that it's 0-based;
- **the end:** the end column, as it is, since the last base of a 1-based, inclusive interval is the end of the 0-based, half-open one;
- **for SAM:** the end from the CIGAR, the start plus the bases that its `M`, `D` and `N` operations cover;
- **for VCF:** the end from the length of REF, or from `INFO/END` where it's past the start.

So a feature gives the same interval in a BED file, which is 0-based, and in a GFF file, which is 1-based:

```@example positions
using HTSLib
using HTSLib: bgzip, tabix_index, records

directory = mktempdir()
bed = joinpath(directory, "feature.bed")
gff = joinpath(directory, "feature.gff")
write(bed, "chr1\t1000\t1500\n")
write(gff, "chr1\tsource\tgene\t1001\t1500\t.\t+\t.\tID=gene1\n")

intervals = map((bed, gff)) do path
    compressed = bgzip(path)
    tabix_index(compressed)
    TabixFile(compressed) do file
        record = only(records(file))
        (HTSLib.is_zero_based(file), record.start, record.stop, HTSLib.position(record))
    end
end
```

`zero_based` of [`tabix_index`](@ref HTSLib.tabix_index) changes only how htslib reads the begin column: `true`, as `tabix -0` does, reads it as 0-based, and `false` as 1-based.

### Starts before the first base

A start before the first base, 0 in a 1-based file, or a negative one in a 0-based file, isn't a position.
As htslib does, it's read as the first base, and htslib warns, on the standard error:

```
[W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
```

The `-0` option is `zero_based = true`: the warning usually means that a 0-based file was indexed as 1-based.

```@example positions
one_based = joinpath(directory, "starts_at_zero.txt")
write(one_based, "chr1\t0\t5\n")
compressed = bgzip(one_based)
tabix_index(compressed; preset = :bed, zero_based = false)
TabixFile(compressed) do file
    record = only(records(file))
    record.start, record.stop
end
```

## Reading the columns yourself

A record's line is the file's, so its columns keep the file's positions.
To read them, check [`HTSLib.is_zero_based`](@ref), and split the line without copying it, with MemoryViews' `split_each`:

```@example positions
using MemoryViews: split_each
using StringViews: StringView

compressed = gff * ".gz"
TabixFile(compressed) do file
    record = only(records(file))
    columns = collect(split_each(HTSLib.line(record), UInt8('\t')))
    begin_column = parse(Int, StringView(columns[4]))
    start = HTSLib.is_zero_based(file) ? begin_column : begin_column - 1
    begin_column, start == record.start
end
```

## Queries

A query's positions are 1-based and inclusive, as a region string's are, whatever the file's positions:
`overlaps(file, "chr1", 1001, 2000)` is the region `chr1:1001-2000`, which htslib reads as the 0-based, half-open interval from 1000 to 2000.
A record overlaps a region if they share a base: `record.start < right` and `record.stop > left - 1`.

## Alignments

SAM counts from 1: POS is the first base that the CIGAR covers, or 0 for a record without a position.
[`HTSLib.end_position`](@ref) is the last base that the CIGAR covers, as htslib's `bam_endpos` computes it:
- a record with `M`, `D`, `N`, `=` or `X` operations ends at `pos` plus the bases they cover, less one;
- an unmapped record, or one without a CIGAR, covers one base, `pos`;
- a record without a position ends at 0, so that [`HTSLib.interval_0based`](@ref) is `(-1, 0)`, as htslib indexes it.

## GenomicFeatures

GenomicFeatures counts from 1, inclusive, so the extension converts a Tabix record's interval, and takes an alignment's as it is; [The GenomicFeatures extension](@ref) describes the rest of its interface:

| Record | `left_position` | `right_position` |
|:--|:--|:--|
| Tabix | `record.start + 1` | `record.stop` |
| alignment | `HTSLib.pos(record)` | `HTSLib.end_position(record)` |
