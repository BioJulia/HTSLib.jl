# Tabix Index Tests
# =================

@testset "tabix_index" begin
    directory = mktempdir()
    compressed = bgzip(write_lines(directory, "index.bed", BED_LINES))
    @test tabix_index(compressed) == compressed * ".tbi"
    @test tabix_index(compressed; csi = true) == compressed * ".csi"
    @test tabix_index(compressed; csi = true, min_shift = 12, index = joinpath(directory, "other.csi")) == joinpath(directory, "other.csi")
    @test all(isfile, (compressed * ".tbi", compressed * ".csi", joinpath(directory, "other.csi")))

    # Columns, positions, the meta character and skipped lines are passed to htslib.
    custom = bgzip(write_lines(directory, "custom.txt", ["x\tchr1\t11\t20", "x\tchr1\t31\t40"]; header = ["skipped line", "%comment"]))
    index = tabix_index(custom; preset = :gff, sequence_column = 2, begin_column = 3, end_column = 4, zero_based = false, meta_char = '%', line_skip = 1)
    TabixFile(custom; index) do file
        @test [(HTSLib.sequence_name(record), record.start, record.stop) for record in records(file)] == [("chr1", 10, 20), ("chr1", 30, 40)]
    end
    @test HTSLib.is_zero_based(TabixIndex(compressed))
    @test !HTSLib.is_zero_based(TabixIndex(custom; index))
    @test HTSLib.is_zero_based(TabixIndex(custom; index = tabix_index(custom; preset = :gff, sequence_column = 2, begin_column = 3, end_column = 4, zero_based = true, meta_char = '%', line_skip = 1, index = joinpath(directory, "custom_0.tbi"))))

    # The BED preset is 0-based, and zero_based = false makes it 1-based, as `tabix -s1 -b2 -e3` does: the start of 11 is the first base, 10 from 0, and the end of 20 the last.
    one_based = bgzip(write_lines(directory, "one_based.bed", ["chr1\t11\t20", "chr1\t31\t40"]))
    tabix_index(one_based; zero_based = false)
    TabixFile(one_based) do file
        @test !HTSLib.is_zero_based(file)
        @test [(record.start, record.stop) for record in records(file)] == [(10, 20), (30, 40)]
        @test [String(record) for record in overlaps(file, "chr1", 20, 31)] == ["chr1\t11\t20", "chr1\t31\t40"]
    end

    # As htslib does, a start before the first base is read as the first base, with a warning, for 1-based and 0-based files alike.
    before_one = bgzip(write_lines(directory, "before_one.bed", ["chr1\t0\t5"]))
    before_zero = bgzip(write_lines(directory, "before_zero.bed", ["chr1\t-1\t5"]))
    for (path, zero_based) in ((before_one, false), (before_zero, true))
        TabixFile(path; index = tabix_index(path; zero_based)) do file
            @test [(record.start, record.stop) for record in records(file)] == [(0, 5)]
        end
    end

    # The preset is detected as `tabix` detects it, from the name, ignoring case, or else from the contents of a VCF file.
    gff = bgzip(write_lines(directory, "genes.gff", ["chr1\tsource\tgene\t11\t20\t.\t+\t.\tID=a"]))
    upper = bgzip(write_lines(directory, "upper.bed", BED_LINES), joinpath(directory, "UPPER.BED.GZ"))
    vcf = bgzip(write_lines(directory, "calls.vcf", ["chr1\t11\t.\tA\tG\t.\t.\t."]; header = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]))
    unnamed_vcf = cp(vcf, joinpath(directory, "calls.txt.gz"))
    columns(path) = (config = HTSLib.config(TabixIndex(path)); (config.preset, config.sc, config.bc, config.ec))
    for (path, preset) in ((compressed, :bed), (upper, :bed), (gff, :gff), (vcf, :vcf), (unnamed_vcf, :vcf))
        tabix_index(path)
        expected = HTSLib.preset_config(preset)
        @test columns(path) == (expected.preset, expected.sc, expected.bc, expected.ec)
    end
    TabixFile(gff) do file
        @test !HTSLib.is_zero_based(file)
        @test [(record.start, record.stop) for record in records(file)] == [(10, 20)]
    end
    TabixFile(unnamed_vcf) do file
        @test [(record.start, record.stop) for record in records(file)] == [(10, 11)]
    end
    # Where `tabix` would use the GFF preset for a file that it cannot detect, a preset must be given.
    unnamed = bgzip(write_lines(directory, "lines.txt", BED_LINES))
    @test_throws "pass `preset`" tabix_index(unnamed)
    # `tabix` detects BAM, CRAM and BCF files by name or contents, and indexes them with htslib's other indexers, rather than with a preset.
    for (name, indexer) in (("lines.bam", "bam_index_build"), ("lines.CRAM", "bam_index_build"), ("lines.bcf", "bcf_index_build"))
        @test_throws indexer tabix_index(cp(unnamed, joinpath(directory, name)))
    end
    # The header of a BAM file without text or sequences, which htslib detects from its magic number.
    bam = joinpath(directory, "empty.alignments")
    write(bam, "BAM\x01", Int32(0), Int32(0))
    @test_throws "is a BAM file" tabix_index(bgzip(bam))
    @test tabix_index(unnamed; preset = :bed) == unnamed * ".tbi"

    @test_throws ErrorException tabix_index(write_lines(directory, "not_compressed.bed", BED_LINES); preset = :bed)
    @test_throws ArgumentError tabix_index(compressed; preset = :unknown)
    @test_throws ArgumentError tabix_index(compressed; min_shift = 12)
    @test_throws ArgumentError tabix_index(compressed; csi = true, min_shift = 0)

    # TabixIndex lifecycle tests
    idx = TabixIndex(compressed)
    @test sequence_names(idx) == ["chr1", "chr2"]
    close(idx)
    @test !isopen(idx)
    @test_throws ArgumentError sequence_names(idx)
    close(idx)
    @test_throws ErrorException TabixIndex(write_lines(directory, "plain.bed", BED_LINES))
end
