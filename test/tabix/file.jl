# Tabix File Tests
# ================

@testset "TabixFile" begin
    directory = mktempdir()
    plain = write_lines(directory, "file.bed", BED_LINES; header = ["# header"])
    compressed = bgzip(plain)
    tbi = tabix_index(compressed)
    csi = tabix_index(compressed; csi = true)

    for index in (nothing, tbi, csi)
        TabixFile(compressed; index) do file
            @test sequence_names(file) == ["chr1", "chr2"]
            # Positions are 1-based and inclusive; the record [10, 20) covers bases 11 to 20.
            @test [String(record) for record in overlaps(file, "chr1:20-31")] == BED_LINES[1:2]
            @test [String(record) for record in query(file, "chr1:20-31")] == BED_LINES[1:2]
            @test isempty(overlaps(file, "chr1:21-30"))
            @test isempty(overlaps(file, "chr1", 1, 10))
            @test [String(record) for record in overlaps(file, "chr1", 1, 11)] == BED_LINES[1:1]
            @test [String(record) for record in overlaps(file, "chr1", 20, 31)] == BED_LINES[1:2]
            @test overlaps(file, "chr1", 20:31) == overlaps(file, "chr1", 20, 31)
            @test overlaps(file, "chr2") == [TabixRecord("chr2", 0, 5, BED_LINES[3])]
            @test isempty(overlaps(file, "chrUnknown"))
            @test isempty(overlaps(file, "chrUnknown", 1, 100))
        end
    end

    TabixFile(compressed; threads = 2) do file
        @test String(only(overlaps(file, "chr2"))) == BED_LINES[3]
        @test occursin("file.bed.gz", repr(file))
    end

    # A plain file is read with the index of its compressed copy, but not queried.
    TabixFile(plain; index = tbi) do file
        @test [String(record) for record in records(file)] == BED_LINES
        @test_throws ArgumentError overlaps(file, "chr1")
        @test_throws ArgumentError overlaps(file, "chr1", 1, 10)
    end

    file = TabixFile(compressed)
    close(file)
    @test !isopen(file)
    @test occursin("closed", repr(file))
    @test_throws ArgumentError overlaps(file, "chr1")
    @test_throws ArgumentError each_record(file)
    close(file)

    @test_throws ErrorException TabixFile(joinpath(directory, "missing.bed.gz"))

    # A record that htslib cannot parse is reported with its line.
    # The index is of another file, which htslib warns of if that file is older.
    invalid = bgzip(write_lines(directory, "invalid.bed", ["chr1\t10\t20", "not a record"]))
    TabixFile(invalid; index = tbi) do file
        @test_throws ErrorException records(file)
    end
end
