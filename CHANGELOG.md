# Changelog

All notable changes to this project are documented in this file.

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- Low-level bindings in `HTSLib.Bindings`: a `ccall` of each C function of htslib that the package uses, and others for packages that need them, such as writing alignments and building their indexes, with the C structs and constants.
- A generator of the low-level bindings, `gen/generate.jl`, which writes their functions, plain structs and constants with Clang.jl from the headers of htslib_jll, as `gen/selection.jl` selects them, and checks the layouts of the structs written by hand.
- htslib_jll's compat is `~1.19.1`, the version the bindings are generated from; the generator fails if the compat admits another minor version.
- `HTSLib.version` and `HTSLib.is_available`.
- Requires Julia 1.12, for the `public` keyword and workspaces.
- The tests are a project of the package's workspace, resolved together with it in one manifest; the generator of the bindings, `gen/`, is a project of its own.
