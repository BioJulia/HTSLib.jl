# Changelog

All notable changes to this project are documented in this file.

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- Low-level bindings in `HTSLib.Bindings`: a `ccall` of each C function of htslib that the package uses, and others for packages that need them, such as writing alignments and building their indexes, with the C structs and constants.
- A generator of the low-level bindings, `gen/generate.jl`, which writes their functions, plain structs and constants with Clang.jl from the headers of htslib_jll, as `gen/selection.jl` selects them, and checks the layouts of the structs written by hand.
- htslib_jll's compat is `~1.19.1`, the version the bindings are generated from; the generator fails if the compat admits another minor version.
- `HTSLib.version` and `HTSLib.is_available`.
- BGZF compression with `bgzip`, and decompression with `bgzip_decompress`.
- Tabix indexes, TBI and CSI, built with `tabix_index`, which detects the preset from the file as `tabix` does, and throws for BAM, CRAM and BCF files, which `tabix` indexes with htslib's other indexers, and loaded with `TabixIndex` and `TabixFile`; `HTSLib.is_zero_based` tells whether an index is for 0-based positions.
- Records read by region with `overlaps`, or its alias `query`, and `each_overlap`, and from start to end with `records` and `each_record`; each query and each read opens the file again, so that they may interleave.
- The types are exported, and the functions are public, imported by name, `using HTSLib: overlaps`, or qualified, since their common names would clash with those of packages loaded alongside.
- Tabix records, `TabixRecord`, whose sequence name `sequence_name` gives, after the tabix specification's "Column for the sequence name", and whose line `line` gives as an `ImmutableMemoryView{UInt8}`.
- SAM, BAM and CRAM files read with `AlignmentFile`, by region with `overlaps`, or its alias `query`, and `each_overlap`, and from start to end with `records` and `each_record`, which open the file again, as Tabix queries do.
- Alignment records, `AlignmentRecord`, with their SAM lines and the names and positions that htslib gives them, through accessors named after the SAM specification's fields, `qname`, `rname` and `pos`, with the specification's other names for them as aliases, `template_name`, `reference_name` and `position`, and `end_position`, the last base of the CIGAR; `qname` is a `StringView` of the record's data, and `rname` the header's `String`, read once, so that neither allocates.
- Records that own their bytes, in Julia's memory, `TabixRecord{HTSLib.Owned}` and `AlignmentRecord{HTSLib.Owned}`, or borrow a buffer of htslib, `TabixRecordView` and `AlignmentRecordView`, which are `TabixRecord{HTSLib.Borrowed}` and `AlignmentRecord{HTSLib.Borrowed}`; all are immutable, and a view equals its copy.
- Iterators that yield owned records by default, and with `copy = false`, views, which allocate nothing; `copy` keeps a view as an owned record, allocating only its bytes, and `collect` copies every view.
- A GenomicFeatures extension, which gives alignments and Tabix records the interval interface, `group_name`, `left_position`, `right_position` and `strand`, and converts them to `GenomicInterval`s, 1-based and inclusive, with the record as the metadata of `GenomicInterval(record)`, or none for `GenomicInterval{Nothing}(record)`.
- Requires Julia 1.12, for the `public` keyword and workspaces.
- The tests are a project of the package's workspace, resolved together with it in one manifest; the generator of the bindings, `gen/`, is a project of its own.
