"""
    HTSLib

Julia bindings to [htslib](https://github.com/samtools/htslib), the C library behind `bgzip`, `tabix` and samtools, from htslib_jll.

HTSLib compresses files with BGZF, builds and loads tabix indexes, and reads records by region or from start to end, as the `bgzip` and `tabix` command-line tools do.
Its use is as a reference implementation, to validate and benchmark Julia packages against htslib in-process, so it reports what htslib finds rather than offering an API of its own for the formats.
The low-level bindings are in [`HTSLib.Bindings`](@ref).
"""
module HTSLib

import htslib_jll
using MemoryViews: ImmutableMemoryView
using StringViews: StringView

# The types are exported, and the functions are public: they are imported by name, `using HTSLib: overlaps`, or qualified, `HTSLib.overlaps`, since their common names would clash with those of packages loaded alongside.
export
    # Abstract types.
    AbstractHTSFile,
    AbstractHTSRecord,

    # Tabix indexes and records.
    TabixIndex,
    TabixFile,
    TabixRecordView,
    TabixRecord

# BGZF compression, tabix indexes, and reading records.
public bgzip, bgzip_decompress, tabix_index, sequence_names, overlaps, query, each_overlap, records, each_record

public Bindings, is_available, version

# Accessors of tabix records.
public position, interval_0based, sequence_name, line

# Whether the positions of a tabix-indexed file are 0-based.
public is_zero_based

# Whether a record owns its bytes or borrows htslib's.
public Owned, Borrowed

"""
    AbstractHTSFile

Abstract supertype for all HTSLib file handles (`TabixFile`).
"""
abstract type AbstractHTSFile end

"""
    AbstractHTSRecord

Abstract supertype for all HTSLib records: `TabixRecord`, which owns its bytes, as a `TabixRecord{HTSLib.Owned}`, or borrows a buffer of htslib, as a `TabixRecordView`, which is `TabixRecord{HTSLib.Borrowed}`.
"""
abstract type AbstractHTSRecord end

include("Bindings.jl")
include("ownership.jl")
include("bgzf.jl")
include("tabix/tabix.jl")

"""
    query(file, region...)

An alias of [`overlaps`](@ref), after htslib's region query, `tbx_itr_querys`, and the region query of `tabix`.
"""
const query = overlaps

"""
    HTSLib.is_available() -> Bool

Whether htslib_jll provides htslib for this platform.
"""
is_available() = htslib_jll.is_available()

"""
    HTSLib.version() -> VersionNumber

The version of the loaded htslib.
"""
version() = VersionNumber(match(r"^\d+(\.\d+){0,2}", Bindings.hts_version()).match)

end # module HTSLib
