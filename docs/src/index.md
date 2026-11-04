# HTSLib.jl

HTSLib.jl calls [htslib](https://github.com/samtools/htslib), the C library behind `bgzip`, `tabix` and samtools, in-process, from Julia.
It loads htslib from [htslib_jll](https://github.com/JuliaBinaryWrappers/htslib_jll.jl), so there is nothing to install besides the package.

Its main use is as a reference implementation: a way to ask htslib what it makes of a file, and to compare a Julia package's answers, or its speed, with htslib's, without starting a process for each call.
So it reports what htslib does, as htslib does it: its positions, its parsing, its corrections and its messages.

## What this package covers

- **BGZF:** compressing and decompressing files, as `bgzip` does.
- **Tabix:** building TBI and CSI indexes, with htslib's presets or columns of your own, as `tabix` does, and reading the records of an indexed file, by region or from start to end.
- **SAM, BAM and CRAM:** reading alignments and their headers, by region with an index, or from start to end, as `samtools view` does.
- **Records:** records that own their bytes, or views of htslib's buffers that read a file without allocating.
- **Low-level bindings:** `HTSLib.Bindings`, a `ccall` of each C function the package uses, and others that packages need, with htslib's structs and constants.

A tabix record is a line with htslib's interval: the package doesn't parse its other columns, whatever the format, BED, GFF, VCF or another.
It doesn't read BCF files, write alignments, or read FASTA files, except as htslib does for CRAM references.

## Installation and platforms

```julia
using Pkg
Pkg.add("HTSLib")
```

The versions of Julia and of htslib that HTSLib.jl supports are the compat bounds of its `Project.toml`.
Its bindings are generated from the one minor version of htslib that those bounds admit, and `HTSLib.version()` gives the version that's loaded.

htslib_jll builds htslib for Linux, macOS and FreeBSD, but not for Windows.
On Windows, the package loads, so that packages that depend on it do too, but `HTSLib.is_available()` is `false`, and every call into htslib fails.

## The manual

- [BGZF and tabix](@ref) compresses files, builds and loads indexes, and reads tabix-indexed files.
- [SAM, BAM and CRAM](@ref) reads alignment files, their headers and their records.
- [Positions and intervals](@ref) describes 0-based and 1-based positions throughout: in files, in records, in queries, and as htslib corrects them.
- [Records and iteration](@ref) describes records that own their bytes and views that borrow htslib's, and how to read without allocating.
- [The GenomicFeatures extension](@ref) gives records the interval interface of GenomicFeatures.
- [Messages and errors](@ref) describes what htslib reports, and what the package throws.
- The [Reference](@ref) and [Low-level bindings](@ref) list every documented function and type.

## First steps

A BED file, compressed and indexed, then read by region:

```@example index
using HTSLib
using HTSLib: bgzip, tabix_index, overlaps

directory = mktempdir()
path = joinpath(directory, "features.bed")
write(path, "chr1\t1000\t1500\tfirst\nchr1\t1999\t2100\tsecond\nchr2\t0\t50\tthird\n")

compressed = bgzip(path)
tabix_index(compressed)

TabixFile(compressed) do file
    overlaps(file, "chr1:1001-2000")
end
```

Each record has its sequence, its interval and its line.
The interval is htslib's, 0-based and half-open, whatever the file's own positions are; [Positions and intervals](@ref) describes how each part of the package counts.

The example's second line imports the functions it calls by name.
HTSLib.jl exports its types, such as `TabixFile` and `AlignmentRecord`, but none of its functions, since their common names would clash with those of packages loaded alongside, such as GenomicFeatures' `each_overlap` and htslib_jll's `bgzip`.
So a program imports the functions it calls, as the example does, or qualifies them, as in `HTSLib.overlaps(file, region)`.
