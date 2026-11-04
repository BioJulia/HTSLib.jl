# BAM Format Tests
# ================

@testset "BAM Format" begin
    directory = mktempdir()
    sam = write_lines(directory, "records.sam", SAM_LINES; header = [chomp(SAM_HEADER)])
    bam = convert_alignments(sam, joinpath(directory, "records.bam"), "wb")

    @testset "read BAM" begin
        AlignmentFile(bam) do file
            @test HTSLib.format(file) == :bam
            @test sequence_names(file) == ["chr1", "chr2"]
            @test sequence_lengths(file) == [1000, 500]
            @test occursin("2 reference sequences", repr(file.header))
            @test names_of(records(file)) == ["read1", "read2", "read3", "placed", "unplaced"]
            @test all(record -> record isa AlignmentRecord{HTSLib.Owned}, each_record(file))
            @test String(file.header) == SAM_HEADER
            @test string.(records(file)) == SAM_LINES
        end
    end

    @testset "query $(index_type) index" for (min_shift, index_type) in [(0, "BAI"), (14, "CSI")]
        index = index_alignments(bam; min_shift = min_shift)
        AlignmentFile(bam; index = index) do file
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

    @testset "lifecycle" begin
        file = AlignmentFile(bam)
        @test occursin("records.bam", repr(file))
        close(file)
        @test !isopen(file)
        @test occursin("closed", repr(file))
        @test_throws ArgumentError overlaps(file, "chr1")
        @test_throws ArgumentError each_record(file)
        @test_throws ArgumentError HTSLib.format(file)
        close(file)

        unindexed = cp(bam, joinpath(directory, "unindexed.bam"))
        AlignmentFile(unindexed) do unindexed_file
            @test_throws ErrorException overlaps(unindexed_file, "chr1")
        end
        # The copy is newer than the index, which htslib warns of.
        AlignmentFile(unindexed; index = bam * ".bai") do indexed_file
            @test names_of(overlaps(indexed_file, "chr1")) == ["read1", "read2"]
        end
        @test_throws ErrorException AlignmentFile(joinpath(directory, "missing.bam"))
    end
end
