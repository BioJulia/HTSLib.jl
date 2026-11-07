# Low-level bindings

`HTSLib.Bindings` has a `ccall` of each C function that HTSLib.jl uses, and the C structs, most of them generated with Clang.jl from the headers of htslib_jll.
Each function carries the C declaration it binds in a comment in the source; the structs written by hand, and the module, are documented here.

```@autodocs
Modules = [HTSLib.Bindings]
Order = [:module, :type, :constant, :function]
```
