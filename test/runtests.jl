# Tests of HTSLib.jl
# ==================
#
# The tests check the bindings and the thin layer over them, not htslib itself.
# The bindings must match the C structs and name symbols that libhts has.
# The layer must pass options and positions to htslib correctly, iterate and free what it allocates, and report errors.
# Each test uses a few records whose expected results are worked out by hand.

using Test
using HTSLib
using HTSLib: bgzip, bgzip_decompress, tabix_index, sequence_names, overlaps, query, each_overlap, records, each_record
using MemoryViews: ImmutableMemoryView
using StringViews: StringView
# The functions of the Libdl standard library, which Base has without the dependency.
using Base.Libc.Libdl: dlsym
import htslib_jll

# Write `lines` to `name` in `directory`, after a header, and return the path.
function write_lines(directory, name, lines; header = String[])
    path = joinpath(directory, name)
    write(path, join(vcat(header, lines), '\n'), '\n')
    return path
end

# The C symbols that Bindings binds, from its `ccall`s, `cglobal`s and `@ccall`s.
function bound_symbols()
    paths = [joinpath(pkgdir(HTSLib), "src", "Bindings.jl"); joinpath.(pkgdir(HTSLib), "src", "Bindings", ["types.jl", "functions.jl"])]
    source = join(read.(paths, String))
    names = [m[1] for m in eachmatch(r"\(:([A-Za-z0-9_]+), libhts\)", source)]
    append!(names, [m[1] for m in eachmatch(r"libhts\.([A-Za-z0-9_]+)\(", source)])
    return unique(names)
end

const BED_LINES = ["chr1\t10\t20\ta", "chr1\t30\t40\tb", "chr2\t0\t5\tc"]

if !HTSLib.is_available()
    @warn "htslib_jll is not available on this platform; skipping the tests."
else

    @testset "HTSLib" begin
        directory = mktempdir()

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

        @testset "bgzip" begin
            plain = write_lines(directory, "bgzip.bed", BED_LINES)
            @test bgzip(plain) == plain * ".gz"
            @test isfile(plain)
            # More than the 1 MiB buffer of `bgzip`, so that it writes several times.
            large = write_lines(directory, "large.txt", [string(k) for k in 1:300_000])
            for (keywords, name) in [((;), "default.gz"), ((; level = 0), "level0.gz"), ((; level = 9), "level9.gz"), ((; threads = 2), "threads.gz")]
                destination = bgzip(large, joinpath(directory, name); keywords...)
                @test read(`$(htslib_jll.bgzip()) -dc $(destination)`, String) == read(large, String)
                decomp = bgzip_decompress(destination, joinpath(directory, name * ".out"); threads = 2)
                @test read(decomp, String) == read(large, String)
            end
            @test filesize(joinpath(directory, "level0.gz")) > filesize(joinpath(directory, "level9.gz"))
            @test_throws ArgumentError bgzip(plain, joinpath(directory, "invalid.gz"); level = 10)
        end

        @testset "Type Hierarchy" begin
            @test TabixFile <: AbstractHTSFile
            # A record owns its bytes or borrows htslib's, which its type parameter says; a view is a borrowing record.
            @test TabixRecord <: AbstractHTSRecord
            @test TabixRecordView === TabixRecord{HTSLib.Borrowed} && TabixRecordView <: TabixRecord
        end

        include("tabix/tabix.jl")
    end

end
