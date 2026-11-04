# HTSLib.jl

Julia bindings to [htslib](https://github.com/samtools/htslib), the C library behind `bgzip`, `tabix` and samtools, loaded from [htslib_jll](https://github.com/JuliaBinaryWrappers/htslib_jll.jl).

Currently, this package's main use is as a reference implementation in-process, to validate and benchmark Julia packages against htslib without starting a process for each call.
Pull requests to fill the gap for other uses are welcome.

## Installation

```julia
using Pkg
Pkg.add("HTSLib")
```

htslib_jll builds htslib for Linux, macOS and FreeBSD, but not for Windows: there the package loads, but `HTSLib.is_available()` is `false`, and every call into htslib fails.
`HTSLib.version()` gives the version of the loaded htslib.

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
