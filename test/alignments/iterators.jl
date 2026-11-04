# Alignment Iterator Tests
# ========================

@testset "Alignment Iterators" begin
    directory = mktempdir()
    sam = write_lines(directory, "records.sam", SAM_LINES; header = [chomp(SAM_HEADER)])
    bam = convert_alignments(sam, joinpath(directory, "records.bam"), "wb")
    index_alignments(bam; min_shift = 0)

    file = AlignmentFile(bam)

    # Iterators yield records of their own, or views with `copy = false`, which `collect` copies.
    @test eltype(each_record(file)) == eltype(each_overlap(file, "chr1")) == AlignmentRecord{HTSLib.Owned}
    @test eltype(each_record(file; copy = false)) == eltype(each_overlap(file, "chr1"; copy = false)) == AlignmentRecordView
    @test all(record -> record isa AlignmentRecordView, each_overlap(file, "chr1:119-200"; copy = false))
    @test collect(each_record(file; copy = false)) == records(file)
    @test collect(each_overlap(file, "chr1:119-200"; copy = false)) == overlaps(file, "chr1:119-200")

    # Reading every record opens the file again, so it may interleave with queries.
    overlapping = each_overlap(file, "chr1", 1, 1000; copy = false)
    sequential = each_record(file; copy = false)
    @test HTSLib.qname(first(iterate(overlapping))) == "read1"
    @test HTSLib.qname(first(iterate(sequential))) == "read1"
    kept = copy(first(iterate(overlapping)))
    @test names_of(collect(copy(record) for record in sequential)) == ["read2", "read3", "placed", "unplaced"]

    # Queries of the same file interleave, each reading through a handle of its own: a query that another runs within reads on from where it was, and its view stays valid until it moves on.
    chr1 = each_overlap(file, "chr1", 1, 1000; copy = false)
    view = first(chr1)
    @test names_of(collect(copy(record) for record in each_overlap(file, "chr2"; copy = false))) == ["read3", "placed"]
    @test HTSLib.qname(view) == "read1"
    @test names_of(collect(copy(record) for record in chr1)) == ["read2"]

    # An unfinished query does not disturb the next.
    for record in each_overlap(file, "chr2")
        break
    end
    @test names_of(overlaps(file, "chr1")) == ["read1", "read2"]

    # Copied records stay valid after their file is closed.
    close(file)
    @test (HTSLib.qname(kept), HTSLib.rname(kept)) == ("read2", "chr1")
end
