# CRAM Format Tests
# =================

@testset "CRAM Format" begin
    directory = mktempdir()
    fasta = write_lines(directory, "reference.fa", [">chr1", "ACGT"^250, ">chr2", "ACGT"^125])
    sam = write_lines(directory, "records.sam", SAM_LINES; header = [chomp(SAM_HEADER)])
    cram = convert_alignments(sam, joinpath(directory, "records.cram"), "wc"; reference = fasta)

    @testset "read CRAM" begin
        AlignmentFile(cram; reference = fasta) do file
            @test HTSLib.format(file) == :cram
            @test sequence_names(file) == ["chr1", "chr2"]
            @test sequence_lengths(file) == [1000, 500]
            @test occursin("2 reference sequences", repr(file.header))
            @test names_of(records(file)) == ["read1", "read2", "read3", "placed", "unplaced"]
            @test all(record -> record isa AlignmentRecord{HTSLib.Owned}, each_record(file))
        end
    end

    @testset "query indexed CRAM" begin
        index = index_alignments(cram; min_shift = 0)
        AlignmentFile(cram; index = index, reference = fasta) do file
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
            # A CRAM index refers to the handle it was loaded for, so each query loads it for its own, and queries interleave as they do in other formats.
            chr1 = each_overlap(file, "chr1")
            @test HTSLib.qname(first(chr1)) == "read1"
            @test names_of(overlaps(file, "chr2")) == ["read3", "placed"]
            @test names_of(collect(chr1)) == ["read2"]
        end
    end

    @testset "error handling" begin
        @test_throws ArgumentError AlignmentFile(fasta)
    end
end
