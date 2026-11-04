# HTSLib.jl

Julia bindings to [htslib](https://github.com/samtools/htslib), the C library behind `bgzip`, `tabix` and samtools, loaded from [htslib_jll](https://github.com/JuliaBinaryWrappers/htslib_jll.jl).

HTSLib.jl compresses files with BGZF, builds and loads tabix indexes, and reads records by region or from start to end, as the `bgzip` and `tabix` command-line tools do.
It reads SAM, BAM and CRAM files, from start to end and by region, as samtools does.

Currently, this package's main use is as a reference implementation in-process, to validate and benchmark Julia packages against htslib without starting a process for each call.
Pull requests to fill the gap for other uses are welcome.

## Installation

```julia
using Pkg
Pkg.add("HTSLib")
```

htslib_jll builds htslib for Linux, macOS and FreeBSD, but not for Windows: there the package loads, but `HTSLib.is_available()` is `false`, and every call into htslib fails.
`HTSLib.version()` gives the version of the loaded htslib.

## BGZF and tabix

```julia
using HTSLib
using HTSLib: bgzip, tabix_index, sequence_names, overlaps, records, each_overlap

compressed = bgzip("features.bed")                 # "features.bed.gz", keeping "features.bed"
index = tabix_index(compressed)                    # "features.bed.gz.tbi", with the BED preset, from the name
tabix_index(compressed; csi = true)                # "features.bed.gz.csi"

TabixFile(compressed) do file
    sequence_names(file)                           # ["chr1", "chr2", …]
    overlaps(file, "chr1:1001-2000")               # Vector{TabixRecord}, as `tabix features.bed.gz chr1:1001-2000`
    overlaps(file, "chr1", 1001, 2000)             # the same region, 1-based and inclusive
    records(file)                                  # every record, from start to end
end
```

HTSLib.jl exports its types, such as `TabixFile`, but none of its functions, since their common names would clash with those of packages loaded alongside, so the examples import the functions they call by name; qualified names, such as `HTSLib.overlaps`, work too.

A `TabixRecord` has the sequence name, the interval, and the line of a record; `HTSLib.sequence_name(record)` gives the sequence name, after the tabix specification's "Column for the sequence name".
As in htslib, intervals are 0-based and half-open, whatever the positions of the file: `start` is the position of the first base, and `stop` the position after the last.
The line keeps the positions of the file, and `HTSLib.is_zero_based(file)` tells whether they are 0-based, as the BED preset and `tabix -0` make them.

`query` is an alias of `overlaps`, after htslib's region queries and the `tabix` command.
`each_overlap` and `each_record` iterate over the records that `overlaps` and `records` collect, and yield them as they do.

A record either owns its line, in Julia's memory, or borrows a buffer of htslib, which its type parameter says: a `TabixRecord{HTSLib.Owned}`, or a `TabixRecordView`, which is `TabixRecord{HTSLib.Borrowed}`.
To read without allocating, `copy = false` yields views, which are only valid until the next iteration; `copy(record)` keeps one, and `collect` copies them all.
`HTSLib.line(record)` gives the line as an `ImmutableMemoryView{UInt8}`, and `StringView(HTSLib.line(record))`, of StringViews, gives it as a string without copying it:

```julia
TabixFile(compressed) do file
    total = 0
    for record in each_overlap(file, "chr1:1001-2000"; copy = false)
        total += record.stop - record.start        # String(record) copies the line
    end
    total
end
```

| Function | `tabix` and `bgzip` equivalent |
|:--|:--|
| `bgzip(source, destination; level, threads)` | `bgzip -c -l level -@ threads source > destination` |
| `tabix_index(path; preset, csi, min_shift, …)` | `tabix -p preset -C -m min_shift path` |
| `TabixFile(path; index, threads)` | the file and index that `tabix path` opens |
| `overlaps(file, region)`, `each_overlap(file, region)` | `tabix path region` |
| `records(file)`, `each_record(file)` | `bgzip -dc path`, without header lines |

`tabix_index` also takes the columns and header settings of `tabix -s`, `-b`, `-e`, `-0`, `-c` and `-S`; see its docstring.
Without `preset`, it detects the preset as `tabix` does, from the name or the contents of the file, and where `tabix` would use the GFF preset for a file that it cannot detect, it throws an `ArgumentError` that asks for one; so do BAM, CRAM and BCF files, which `tabix` indexes with htslib's other indexers.

## SAM, BAM and CRAM

```julia
using HTSLib
using HTSLib: sequence_names, sequence_lengths, records, overlaps

AlignmentFile("reads.bam") do file
    HTSLib.format(file)                            # :bam
    sequence_names(file), sequence_lengths(file)   # the reference sequences of the header
    String(file.header)                            # the header text
    records(file)                                  # every record, as `samtools view reads.bam`
    for record in overlaps(file, "chr1:1001-2000") # as `samtools view reads.bam chr1:1001-2000`
        string(record)                             # the SAM line, as samtools prints it
        HTSLib.pos(record), HTSLib.end_position(record)
    end
end
```

HTSLib.jl only reads alignment files, to report what htslib finds in them.
Test files can be converted from SAM and indexed with htslib through `HTSLib.Bindings`, as samtools does, which the tests of HTSLib.jl show.

A CRAM file needs its FASTA reference, `AlignmentFile("reads.cram"; reference = "genome.fa")`, or htslib looks it up as samtools does, which may download it.
As samtools does, htslib adds `MD` and `NM` tags computed from the reference to the records of a CRAM file.

As Tabix records are, alignment records are owned, `AlignmentRecord{HTSLib.Owned}`, or views, `AlignmentRecordView`, which `copy = false` yields.
`HTSLib.data(record)` gives the variable-length fields as htslib lays them out, as an `ImmutableMemoryView{UInt8}`.
`string(record)` gives the SAM line of a record, with every field.
The accessors give what htslib makes of a record's name and position.
They take the names of the SAM specification's fields, with aliases where the specification names a field otherwise, and, like the functions, are public but not exported:

| Accessor | Value | Aliases |
|:--|:--|:--|
| `HTSLib.qname` | QNAME, "Query template NAME", as a `StringView` of the record's data, valid as long as the record is | `template_name` |
| `HTSLib.rname` | RNAME, "Reference sequence NAME", as the header's `String`, or `nothing` without one | `reference_name` |
| `HTSLib.pos` | POS, "1-based leftmost mapping POSition", or 0 without one | `position` |
| `HTSLib.end_position` | the last reference base that the CIGAR covers, 1-based, as htslib's `bam_endpos` computes it: an unmapped record, or one without a CIGAR, covers one base, and a record without a position ends at 0 | — |
| `HTSLib.reference_index` | the 1-based index of RNAME into `sequence_names`, or 0 without one | — |

| Function | samtools equivalent |
|:--|:--|
| `AlignmentFile(path; index, reference)`, `records(file)`, `each_record(file)` | `samtools view path` |
| `overlaps(file, region)`, `each_overlap(file, region)` | `samtools view path region`, with `"*"` for records without a position |

## GenomicFeatures

With GenomicFeatures loaded, alignments and Tabix records have its interval interface: `group_name` is `rname` of an alignment and `sequence_name` of a Tabix record, `left_position` and `right_position` are 1-based and inclusive, and `strand` is an alignment's, from its FLAG, or unknown.
`GenomicInterval(record)` converts a record to a `GenomicInterval`, with the record as its metadata, and `GenomicInterval{Nothing}(record)` without.

## Low-level bindings

`HTSLib.Bindings` has a `ccall` of each C function that HTSLib.jl uses, and of others that packages need, such as writing alignments and building their indexes, with htslib's structs and constants.
Each function carries the C declaration it binds in a comment, and structs used only through pointers, such as `BGZF`, `hts_idx_t` and `sam_hdr_t`, are named opaque types.

Most of the bindings are generated with [Clang.jl](https://github.com/JuliaInterop/Clang.jl) from the headers of htslib_jll, as `gen/selection.jl` selects them, and the generator checks the structs written by hand against htslib; [`gen/README.md`](gen/README.md) describes how to add a function and regenerate them.

## Testing

The tests check the bindings and the thin layer over them, not htslib itself: the layouts of the structs, that libhts has every bound symbol, and that the layer passes options and positions to htslib, iterates, frees what it allocates and reports errors, on a few records with results worked out by hand:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Alternatives and inspiration

Alternatives and inspiration for this package include [htslib_jll](https://github.com/JuliaBinaryWrappers/htslib_jll.jl) and [https://github.com/bicycle1885/htslib.jl](https://github.com/bicycle1885/htslib.jl).
