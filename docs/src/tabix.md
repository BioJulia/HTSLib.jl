# BGZF and tabix

A tabix-indexed file is a text file of one record per line, sorted by sequence and start, compressed with BGZF, with an index of the positions of its records.
BGZF is gzip in independent blocks, so `gunzip` reads it, and a reader can seek to a block and decompress from there; the index maps each region to the blocks that hold its records.
This page goes from a plain text file to queries of it, as `bgzip` and `tabix` do.

```@example tabix
using HTSLib
using HTSLib: bgzip, tabix_index, sequence_names, overlaps, each_overlap, records

directory = mktempdir()
lines = ["chr1\t1000\t1500\tfirst", "chr1\t1999\t2100\tsecond", "chr1\t5000\t5001\tthird", "chr2\t0\t50\tfourth"]
path = joinpath(directory, "features.bed")
write(path, "# a header line\n", join(lines, '\n'), '\n')
nothing # hide
```

## Compression

[`bgzip`](@ref HTSLib.bgzip) compresses a file with BGZF, as `bgzip -c` does, and keeps the original:

```@example tabix
compressed = bgzip(path)
basename(compressed)
```

The destination is `source * ".gz"` unless given, `level` is the compression level, from 0, for none, to 9, and `threads` is the number of threads that compress.
[`bgzip_decompress`](@ref HTSLib.bgzip_decompress) does the reverse, as `bgzip -dc` does, to the name without `.gz`, or to `source * ".decompressed"`.

Tabix needs BGZF rather than gzip: htslib can't index a file compressed with `gzip`, since it can't seek within it.

## Building an index

[`tabix_index`](@ref HTSLib.tabix_index) builds an index of a compressed, sorted file, and returns its path:

```@example tabix
basename(tabix_index(compressed))
```

### Presets

An index records which columns hold the sequence name and the interval of each record, and whether the positions are 0-based.
htslib has a preset for each format it knows, as `tabix -p` does:

| `preset` | Sequence | Begin | End | Positions | Comment lines |
|:--|:--|:--|:--|:--|:--|
| `:bed` | 1 | 2 | 3 | 0-based, half-open | `#` |
| `:gff` | 1 | 4 | 5 | 1-based, inclusive | `#` |
| `:vcf` | 1 | 2 | from the length of REF, or `INFO/END` | 1-based | `#` |
| `:sam` | 3 | 4 | from the CIGAR | 1-based | `@` |
| `:psltbl` | 15 | 17 | 18 | 0-based, half-open | `#` |

Without `preset`, `tabix_index` picks it as `tabix` does:
1. from the name, ignoring case: `.bed.gz`, `.gff.gz`, `.vcf.gz` or `.sam.gz`;
2. then from the contents, if htslib detects a VCF file.

Where `tabix` would fall back to GFF for a file it can't detect, `tabix_index` throws an `ArgumentError` that asks for a preset, so that a file isn't indexed with the wrong columns without notice.
BAM, CRAM and BCF files, which `tabix` indexes with htslib's other indexers, `bam_index_build` and `bcf_index_build`, throw an `ArgumentError` too.

### Options

The keywords choose the preset, replace its columns and settings, and choose the kind of index, as `tabix`'s options do:

| Keyword | `tabix` option | Meaning |
|:--|:--|:--|
| `preset` | `-p` | the preset to start from |
| `sequence_column` | `-s` | the 1-based column of the sequence name |
| `begin_column` | `-b` | the 1-based column of the start |
| `end_column` | `-e` | the 1-based column of the end |
| `zero_based` | `-0` | `true` for 0-based, half-open positions; `false` for 1-based, inclusive ones |
| `meta_char` | `-c` | the first character of header and comment lines |
| `line_skip` | `-S` | the number of lines to skip at the start |
| `csi` | `-C` | build a CSI index instead of a TBI index |
| `min_shift` | `-m` | the size of the smallest bins of a CSI index, `2^min_shift` bases |
| `index` | — | the path of the index, instead of the file's path with `.tbi` or `.csi` |
| `threads` | `-@` | the number of threads that decompress the file |

A file of 1-based positions in columns 2 to 4, after a sample column, with `%` before its comments:

```@example tabix
custom = joinpath(directory, "custom.txt")
write(custom, "%sample\tsequence\tbegin\tend\n", "a\tchr1\t1001\t1500\n", "b\tchr1\t2000\t2100\n")
custom = bgzip(custom)
tabix_index(custom; preset = :gff, sequence_column = 2, begin_column = 3, end_column = 4, meta_char = '%')
TabixFile(records, custom)
```

The GFF preset is 1-based, so the start of 1001 is htslib's 1000; [Positions and intervals](@ref) describes how htslib reads each.

### TBI and CSI

A TBI index, the default, holds positions up to 2^29, about 537 million bases; htslib refuses to store a record past that in one, and says to use a CSI index.
`csi = true` builds a CSI index, whose bins at the finest level are `2^min_shift` bases, 16384 by default, and which holds longer sequences.

## Opening a file

[`TabixFile`](@ref) opens a compressed file and loads its index, found next to the file as `tabix` finds it, a `.tbi` or `.csi` index, unless `index` gives its path.
[`TabixIndex`](@ref) loads an index on its own, without its file.

```@example tabix
file = TabixFile(compressed)
sequence_names(file), HTSLib.is_zero_based(file)
```

Close a file with `close`, or pass a function, `TabixFile(f, path)`, which is called with the file, and closes it afterwards; the garbage collector closes a file that's left open.

A plain, uncompressed file can be opened with the index of its compressed copy, and read from start to end, but not queried, since a query seeks to the blocks of a BGZF file.

## Reading from start to end

[`records`](@ref HTSLib.records) reads every record, and [`each_record`](@ref HTSLib.each_record) iterates over them:

```@example tabix
records(file)
```

As when the index was built, the first `line_skip` lines, and lines that start with the comment character, are skipped, so the file's header line, `# a header line`, isn't a record.
A line that htslib can't parse, or whose sequence the index doesn't have, throws an error with its line number; for a line it can't parse, htslib prints its own message first.

## Queries

[`overlaps`](@ref HTSLib.overlaps) collects the records that overlap a region, and [`each_overlap`](@ref HTSLib.each_overlap) iterates over them; `query` is an alias of `overlaps`, after htslib's region queries.
A region is a string, as `tabix file region` takes it, or a sequence name with 1-based, inclusive positions:

```@example tabix
overlaps(file, "chr1", 2000, 2100) == overlaps(file, "chr1:2000-2100") == overlaps(file, "chr1", 2000:2100)
```

htslib parses region strings, so they have its forms:

| Region | Records that overlap |
|:--|:--|
| `"chr1"` | the whole of `chr1` |
| `"chr1:1001-2000"` | bases 1001 to 2000 of `chr1`, 1-based and inclusive |
| `"chr1:1,001-2,000"` | the same; commas in numbers are ignored |
| `"chr1:2000"` | from base 2000 to the end of `chr1` |
| `"chr1:-1000"` | from the start of `chr1` to base 1000 |
| `"{chr1:x}"` | the whole of a sequence whose name has a colon, `chr1:x` |
| `"."` | the whole file |

```@example tabix
[String(record) for record in overlaps(file, "chr1:2000")]
```

A record overlaps a region if it shares a base with it: for the region from `left` to `right`, 1-based and inclusive, a record overlaps if `record.start < right` and `record.stop > left - 1`.
So the BED record `chr1 1999 2100`, whose first base is 2000, overlaps `chr1:2000-2000`, but not `chr1:1999-1999`:

```@example tabix
length(overlaps(file, "chr1:2000-2000")), length(overlaps(file, "chr1:1999-1999"))
```

As in htslib, a sequence that the index doesn't have has no records, and neither does a region that htslib can't parse, after its message.
With a sequence and positions, a `left` below 1 is read as 1, and a `right` below `left` has no records.

Each query, and each read from start to end, opens the file again and reads through a handle of its own, so they may interleave.
Here a query of `chr1` reads its first record, a query of `chr2` runs, and the first reads on with the other two records of `chr1`:

```@example tabix
chr1 = each_overlap(file, "chr1")
first_of_chr1 = String(first(chr1))
of_chr2 = [String(record) for record in overlaps(file, "chr2")]
rest_of_chr1 = [String(record) for record in chr1]
first_of_chr1, of_chr2, rest_of_chr1
```

```@example tabix
close(file)
```

## Equivalents

| HTSLib.jl | `bgzip` and `tabix` |
|:--|:--|
| `bgzip(source; level, threads)` | `bgzip -c -l level -@ threads source > source.gz` |
| `bgzip_decompress(source; threads)` | `bgzip -dc -@ threads source > destination` |
| `tabix_index(path; preset, csi, min_shift)` | `tabix -p preset -C -m min_shift path` |
| `tabix_index(path; sequence_column, begin_column, end_column, zero_based)` | `tabix -s … -b … -e … -0 path` |
| `overlaps(file, region)` | `tabix path region` |
| `records(file)` | `bgzip -dc path`, without header and skipped lines |
| `sequence_names(file)` | `tabix -l path` |
