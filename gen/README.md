# Generating the bindings

The functions, plain structs and constants of `HTSLib.Bindings` are generated with [Clang.jl](https://github.com/JuliaInterop/Clang.jl) from the headers of htslib_jll.
Clang.jl generates bindings for all of htslib's API, more than 350 functions, with htslib's own names and C types.
`selection.jl` selects what `HTSLib.Bindings` binds from that, and sets how its names and types differ.

| File | Contents |
|:--|:--|
| `selection.jl` | The selection: the functions, structs and constants to bind, the rules and overrides for their types and names, and the hand-written structs to check. |
| `generate.jl` | The generator, which writes `src/Bindings/types.jl` and `src/Bindings/functions.jl`. |
| `generator.toml` | Clang.jl's options. |

`src/Bindings.jl` is the module.
It includes the two generated files between the parts written by hand, which are what Clang.jl can't generate, and its opening comment lists the six parts in order.

## Running

```bash
julia --project=gen -e 'using Pkg; Pkg.instantiate()'
```

```bash
julia --project=gen gen/generate.jl
```

The generator takes about 20 seconds.
It then loads `src/Bindings.jl` and checks the offsets and sizes of the structs written there by hand against those that libclang computes, failing if they differ.
With `--full`, it also writes Clang.jl's complete output to `gen/Generated.jl`, which isn't committed, for browsing the rest of htslib's API.

Run it after changing `selection.jl`, and review the diff.

## Updating htslib

The bindings are generated from one minor version of htslib, and HTSLib.jl's compat admits only that version, with a tilde bound, `htslib_jll = "~X.Y.Z"`.
The generator fails if the htslib_jll it runs on isn't within that compat, or if the compat admits another minor version.
To move to a new htslib:
1. set the bound of htslib_jll in both `Project.toml` and `gen/Project.toml` to the new version, and update `gen/`'s manifest;
2. run the generator;
3. review the diff of `src/Bindings/` and the layout check;
4. run the tests of HTSLib.jl and of the packages that use it.

## Adding a function

Add the C function's name to a section of `FUNCTIONS` in `selection.jl`, and regenerate.
The rules at the top of `selection.jl` rewrite Clang.jl's types:
- strings, `const char *`, become `Cstring`;
- out-parameters of `kstring_t`, `int` and `hts_pos_t` become `Ref`s;
- enums become `Cint`;
- the argument name `end` becomes `end_`.

Override what the rules get wrong in the entry, by the C names of arguments:

```julia
:tbx_index_load3 => (arguments = (fnidx = :(Ptr{UInt8}),),),  # an optional string, which takes C_NULL
```

A name that htslib doesn't declare, or an override of an argument that the function doesn't have, is an error.

What Clang.jl can't generate stays in `src/Bindings.jl`, written by hand:
- function-like macros, such as `tbx_itr_queryi`;
- inline functions, such as `sam_itr_next`;
- variadic functions, such as `hts_set_opt`;
- extern globals, such as `tbx_conf_bed`;
- structs with bit fields or unions, such as `bam1_t`, which `LAYOUTS` in `selection.jl` lists for checking.

## Clang.jl's options

- The headers are parsed in the order `bgzf.h`, `hts.h`, `sam.h`, `tbx.h`. `bgzf.h` defines `struct BGZF`, which `hts.h` only declares, and Clang.jl's sanity check fails if the declaration comes first.
- `extract_c_comment_style = "disable"`: extracting doxygen comments makes a run take about 20 minutes, rather than 20 seconds, and the bindings carry their own comments.
