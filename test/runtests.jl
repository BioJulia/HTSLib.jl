# Tests of HTSLib.jl
# ==================
#
# The tests check the bindings and the thin layer over them, not htslib itself.
# The bindings must match the C structs and name symbols that libhts has.
# The layer must pass options and positions to htslib correctly, iterate and free what it allocates, and report errors.
# Each test uses a few records whose expected results are worked out by hand.

using Test
using HTSLib
# The functions of the Libdl standard library, which Base has without the dependency.
using Base.Libc.Libdl: dlsym
import htslib_jll

# The C symbols that Bindings binds, from its `ccall`s, `cglobal`s and `@ccall`s.
function bound_symbols()
    paths = [joinpath(pkgdir(HTSLib), "src", "Bindings.jl"); joinpath.(pkgdir(HTSLib), "src", "Bindings", ["types.jl", "functions.jl"])]
    source = join(read.(paths, String))
    names = [m[1] for m in eachmatch(r"\(:([A-Za-z0-9_]+), libhts\)", source)]
    append!(names, [m[1] for m in eachmatch(r"libhts\.([A-Za-z0-9_]+)\(", source)])
    return unique(names)
end

if !HTSLib.is_available()
    @warn "htslib_jll is not available on this platform; skipping the tests."
else

    @testset "HTSLib" begin
        @testset "Bindings" begin
            # The structs that are read from or passed to htslib have the sizes and offsets of C's on 64-bit platforms.
            if Sys.WORD_SIZE == 64
                Bindings = HTSLib.Bindings
                @test sizeof(Bindings.kstring_t) == 24
                @test sizeof(Bindings.tbx_conf_t) == 24
                @test sizeof(Bindings.tbx_t) == 40
                @test sizeof(Bindings.bam1_core_t) == 48
                @test sizeof(Bindings.bam1_t) == 80
                @test sizeof(Bindings.htsFormat) == 32
                # hts_itr_t and htsFile are longer in C, but their leading fields end with curr_end and the union of streams.
                @test fieldoffset(Bindings.hts_itr_t, 9) == 48
                @test fieldoffset(Bindings.hts_itr_t, 13) == 72
                @test fieldoffset(Bindings.htsFile, 6) == 56
            end
            # Every bound function and global is in libhts.
            symbols = bound_symbols()
            @test length(symbols) > 50
            for symbol in symbols
                @test (symbol => dlsym(htslib_jll.libhts_handle, symbol; throw_error = false) !== nothing) == (symbol => true)
            end
            # A preset global is read as its struct: BED has columns 1, 2 and 3, and 0-based positions.
            bed = HTSLib.Bindings.tbx_conf_bed()
            @test (bed.sc, bed.bc, bed.ec, Char(bed.meta_char)) == (1, 2, 3, '#')
            @test bed.preset & HTSLib.Bindings.TBX_UCSC != 0
            @test HTSLib.version() == VersionNumber(pkgversion(htslib_jll).major, pkgversion(htslib_jll).minor, pkgversion(htslib_jll).patch)
        end
    end

end
