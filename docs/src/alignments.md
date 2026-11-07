# SAM, BAM and CRAM

[`AlignmentFile`](@ref) reads alignment files as samtools does: SAM, plain or compressed with `bgzip`, BAM and CRAM, from start to end, or by region with an index.
htslib detects the format from the file's contents, not its name.

## Making files to read

HTSLib.jl only reads alignment files, but its low-level bindings can write them, as `samtools view` and `samtools index` do.
This page's examples convert a few SAM records to BAM and CRAM files this way:

```@example alignments
using HTSLib
using HTSLib: sequence_names, sequence_lengths, records, each_record, overlaps
const Bindings = HTSLib.Bindings

directory = mktempdir()
sam = joinpath(directory, "reads.sam")
write(sam, """
    @HD\tVN:1.6\tSO:coordinate
    @SQ\tSN:chr1\tLN:1000
    @SQ\tSN:chr2\tLN:500
    read1\t0\tchr1\t100\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII
    read2\t16\tchr1\t200\t30\t10M5D10M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*
    read3\t0\tchr2\t50\t60\t5S20M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTA\t*
    placed\t4\tchr2\t60\t0\t*\t*\t0\t0\tACGT\t*
    unplaced\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*
    """)

# Convert `source` as `samtools view` does: to BAM with mode "wb", or to CRAM with "wc" and a reference.
function convert_alignments(source, destination, mode; reference = nothing)
    input = Bindings.hts_open(source, "r")
    header = Bindings.sam_hdr_read(input)
    output = Bindings.hts_open(destination, mode)
    reference === nothing || Bindings.hts_set_fai_filename(output, reference)
    Bindings.sam_hdr_write(output, header)
    record = Bindings.bam_init1()
    while Bindings.sam_read1(input, header, record) >= 0
        Bindings.sam_write1(output, header, record)
    end
    Bindings.bam_destroy1(record)
    Bindings.sam_hdr_destroy(header)
    Bindings.hts_close(output)
    Bindings.hts_close(input)
    return destination
end

bam = convert_alignments(sam, joinpath(directory, "reads.bam"), "wb")
Bindings.sam_index_build3(bam, C_NULL, 0, 0)    # reads.bam.bai, as `samtools index` builds it
nothing # hide
```

A real program checks what each binding returns; [Low-level bindings](@ref) lists them.

## Opening a file

```@example alignments
file = AlignmentFile(bam)
HTSLib.format(file), sequence_names(file), sequence_lengths(file)
```

`file.header` is the file's [`AlignmentHeader`](@ref), and `String(file.header)` its text:

```@example alignments
print(String(file.header))
```

Close a file with `close`, or pass a function, `AlignmentFile(f, path)`, which is called with the file, and closes it afterwards.
Owned records stay valid after the file is closed, since they keep its header; views don't, as [Records and iteration](@ref) describes.

### CRAM references

A CRAM file stores bases as differences from a reference, so reading it needs the reference's FASTA file, with its `.fai` index next to it, which htslib builds if it's missing:

```@example alignments
fasta = joinpath(directory, "reference.fa")
write(fasta, ">chr1\n", "ACGT"^250, "\n>chr2\n", "ACGT"^125, "\n")
cram = convert_alignments(sam, joinpath(directory, "reads.cram"), "wc"; reference = fasta)
AlignmentFile(cram; reference = fasta) do cram_file
    HTSLib.format(cram_file), length(records(cram_file))
end
```

Without `reference`, htslib looks the reference up as samtools does, through the `REF_PATH` and `REF_CACHE` environment variables, which may download it.
As samtools does, htslib adds `MD` and `NM` tags, computed from the reference, to the records of a CRAM file.

## Reading from start to end

[`records`](@ref HTSLib.records) reads every record, as `samtools view` does, and [`each_record`](@ref HTSLib.each_record) iterates over them:

```@example alignments
[String(HTSLib.qname(record)) for record in each_record(file)]
```

## Queries

A query needs the file's index, a BAI or CSI index of a BAM file or of a SAM file compressed with `bgzip`, or a CRAI index of a CRAM file.
It's found next to the file, as samtools finds it, unless `index` gives its path, and it's loaded at the first query.
A plain SAM file can be read from start to end, but not queried.

[`overlaps`](@ref HTSLib.overlaps) collects the records that overlap a region, as `samtools view file region` does, and [`each_overlap`](@ref HTSLib.each_overlap) iterates over them; `query` is an alias of `overlaps`.

```@example alignments
[String(HTSLib.qname(record)) for record in overlaps(file, "chr1:110-210")]
```

Regions take the forms that tabix queries do, described in [BGZF and tabix](@ref), with sequence names from the header, and two more:

| Region | Records |
|:--|:--|
| `"*"` | the records without a position, at the end of the file |
| `"."` | every record |

```@example alignments
[String(HTSLib.qname(record)) for record in overlaps(file, "*")], length(overlaps(file, "."))
```

A record overlaps a region from `left` to `right` if `HTSLib.pos(record) <= right` and `HTSLib.end_position(record) >= left`.
An unmapped record that's placed, with a position but no CIGAR, as the mate of a mapped read often is, covers one base at its position, as htslib indexes it, so a query of that base finds it, with `read3`, which covers bases 50 to 69:

```@example alignments
[String(HTSLib.qname(record)) for record in overlaps(file, "chr2:60-60")]
```

Each query opens the file again and reads through a handle of its own, so queries may interleave, with each other and with reads from start to end.
A CRAM index refers to the handle it was loaded for, so each query of a CRAM file loads the index again.

## Records

An [`AlignmentRecord`](@ref) holds the fixed-length fields of htslib's `bam1_t`, its variable-length data, and the header of its file.
`string(record)` is its SAM line, as samtools prints it:

```@example alignments
read2 = overlaps(file, "chr1:200-200")[1]
string(read2)
```

The accessors give what htslib makes of the record's name and position, and take the names of the SAM specification's fields, with its other names for them as aliases:

| Accessor | Value | Alias |
|:--|:--|:--|
| [`HTSLib.qname`](@ref) | QNAME, as a `StringView` of the record's data | `template_name` |
| [`HTSLib.rname`](@ref) | RNAME, as the header's `String`, or `nothing` | `reference_name` |
| [`HTSLib.pos`](@ref) | POS, 1-based, or 0 without a position | `position` |
| [`HTSLib.end_position`](@ref) | the last base that the CIGAR covers, 1-based | — |
| [`HTSLib.reference_index`](@ref) | the 1-based index of RNAME in `sequence_names`, or 0 | — |
| [`HTSLib.interval_0based`](@ref) | `(pos - 1, end_position)`, 0-based and half-open | — |

```@example alignments
HTSLib.qname(read2), HTSLib.rname(read2), HTSLib.pos(read2), HTSLib.end_position(read2)
```

`read2`'s CIGAR, `10M5D10M`, covers 25 bases of the reference from 200, with its deletion, so it ends at 224.

The other fields, FLAG, MAPQ, CIGAR, SEQ, QUAL and the tags, are in the SAM line, and in `HTSLib.data(record)`, the variable-length data as htslib lays it out, as bytes.

```@example alignments
close(file)
```

## Equivalents

| HTSLib.jl | samtools |
|:--|:--|
| `records(file)`, `each_record(file)` | `samtools view path` |
| `overlaps(file, region)`, `each_overlap(file, region)` | `samtools view path region` |
| `String(file.header)` | `samtools view -H path` |
| `AlignmentFile(path; reference)` | `samtools view -T reference path` |
| `AlignmentFile(path; index)` | `samtools view -X path index` |
