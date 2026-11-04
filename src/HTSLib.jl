"""
    HTSLib

Julia bindings to [htslib](https://github.com/samtools/htslib), the C library behind `bgzip`, `tabix` and samtools, from htslib_jll.

Its use is as a reference implementation, to validate and benchmark Julia packages against htslib in-process, so it reports what htslib finds rather than offering an API of its own for the formats.
The low-level bindings are in [`HTSLib.Bindings`](@ref).
"""
module HTSLib

import htslib_jll

public Bindings, is_available, version

include("Bindings.jl")

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
