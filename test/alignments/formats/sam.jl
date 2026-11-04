# SAM Format Tests
# ================

@testset "SAM Format" begin
    directory = mktempdir()
    sam = write_lines(directory, "records.sam", SAM_LINES; header = [chomp(SAM_HEADER)])
    sam_gz = bgzip(sam)

    @testset "read plain SAM" begin
        AlignmentFile(sam) do file
            @test HTSLib.format(file) == :sam
            @test sequence_names(file) == ["chr1", "chr2"]
            @test sequence_lengths(file) == [1000, 500]
            @test occursin("2 reference sequences", repr(file.header))
            @test names_of(records(file)) == ["read1", "read2", "read3", "placed", "unplaced"]
            @test all(record -> record isa AlignmentRecord{HTSLib.Owned}, each_record(file))
            @test String(file.header) == SAM_HEADER
            @test string.(records(file)) == SAM_LINES
        end
    end

    @testset "read compressed SAM (.sam.gz)" begin
        AlignmentFile(sam_gz) do file
            @test HTSLib.format(file) == :sam
            @test sequence_names(file) == ["chr1", "chr2"]
            @test sequence_lengths(file) == [1000, 500]
            @test occursin("2 reference sequences", repr(file.header))
            @test names_of(records(file)) == ["read1", "read2", "read3", "placed", "unplaced"]
            @test all(record -> record isa AlignmentRecord{HTSLib.Owned}, each_record(file))
            @test String(file.header) == SAM_HEADER
            @test string.(records(file)) == SAM_LINES
        end
    end

    @testset "query on plain SAM fails" begin
        AlignmentFile(sam) do file
            @test_throws ArgumentError overlaps(file, "chr1")
        end
    end

    @testset "query on indexed compressed SAM (.sam.gz)" begin
        index = index_alignments(sam_gz; min_shift = 14)
        AlignmentFile(sam_gz; index = index) do file
            @test names_of(overlaps(file, "chr1", 119, 199)) == ["read1"]
            @test names_of(query(file, "chr1", 119, 199)) == ["read1"]
            @test isempty(overlaps(file, "chr1", 120, 199))
            @test names_of(overlaps(file, "chr1", 120, 200)) == ["read2"]
            @test names_of(overlaps(file, "chr1:119-200")) == ["read1", "read2"]
            @test names_of(overlaps(file, "chr1", 119:200)) == ["read1", "read2"]
            @test names_of(overlaps(file, "chr2")) == ["read3", "placed"]
            @test names_of(overlaps(file, "*")) == ["unplaced"]
            @test isempty(overlaps(file, "chrUnknown"))
            @test isempty(overlaps(file, "chrUnknown", 1, 100))
        end
    end
end
