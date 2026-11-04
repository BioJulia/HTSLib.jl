# Alignment Record and View Tests
# ===============================

@testset "AlignmentRecord & AlignmentRecordView" begin
    directory = mktempdir()
    sam = write_lines(directory, "records.sam", SAM_LINES; header = [chomp(SAM_HEADER)])
    bam = convert_alignments(sam, joinpath(directory, "records.bam"), "wb")

    (read1, read2, read3, placed, unplaced) = AlignmentFile(records, bam)

    # qname, rname and reference_index
    @test (HTSLib.qname(read1), HTSLib.rname(read1), HTSLib.reference_index(read1)) == ("read1", "chr1", 1)
    # The aliases: the SAM specification's own names for the fields.
    @test (HTSLib.template_name(read1), HTSLib.reference_name(read1), HTSLib.position(read1)) == ("read1", "chr1", 100)
    @test (HTSLib.pos(read1), HTSLib.end_position(read1)) == (100, 119)
    @test HTSLib.interval_0based(read1) == (99, 119)

    @test (HTSLib.pos(read2), HTSLib.end_position(read2)) == (200, 224)
    @test (HTSLib.rname(read3), HTSLib.reference_index(read3), HTSLib.pos(read3), HTSLib.end_position(read3)) == ("chr2", 2, 50, 69)
    # As htslib's bam_endpos computes it, an unmapped record covers one base, and a record without a position ends at 0, as htslib indexes it.
    @test (HTSLib.rname(placed), HTSLib.pos(placed), HTSLib.end_position(placed), HTSLib.interval_0based(placed)) == ("chr2", 60, 60, (59, 60))
    @test (HTSLib.rname(unplaced), HTSLib.reference_index(unplaced), HTSLib.pos(unplaced)) == (nothing, 0, 0)
    @test (HTSLib.end_position(unplaced), HTSLib.interval_0based(unplaced)) == (0, (-1, 0))
    @test HTSLib.reference_name(unplaced) === nothing

    # string formatting and show, of records copied from a file that is now closed
    @test string(read1) == SAM_LINES[1]
    @test occursin("AlignmentRecord(\"read1", repr(read1))

    # AlignmentRecordView accessors
    AlignmentFile(bam) do file
        view1 = first(each_record(file; copy = false))
        @test view1 isa AlignmentRecordView
        @test HTSLib.qname(view1) == "read1"
        @test HTSLib.rname(view1) == "chr1"
        @test HTSLib.reference_index(view1) == 1
        @test HTSLib.reference_name(view1) == "chr1"
        @test HTSLib.pos(view1) == 100
        @test HTSLib.end_position(view1) == 119
        @test HTSLib.position(view1) == 100
        @test HTSLib.interval_0based(view1) == (99, 119)
        @test string(view1) == SAM_LINES[1]
        @test occursin("AlignmentRecordView(\"read1", repr(view1))

        # `copy` gives an immutable record in Julia's memory, whose copy is itself.
        copy1 = copy(view1)
        @test copy1 isa AlignmentRecord{HTSLib.Owned} && !ismutable(copy1)
        @test copy(copy1) === copy1
        @test string(copy1) == string(view1) == SAM_LINES[1]
        @test (HTSLib.qname(copy1), HTSLib.pos(copy1), HTSLib.end_position(copy1)) == ("read1", 100, 119)
        # Copies of one record are equal, with equal hashes; copies of others are not.
        again = copy(first(each_record(file)))
        @test again == copy1 && hash(again) == hash(copy1)
        views = each_record(file; copy = false)
        (_, state) = iterate(views)
        (view2, _) = iterate(views, state)
        second = copy(view2)
        @test HTSLib.qname(second) == "read2" && second != copy1

        # A view equals its copy, and `data` gives its bytes, which can't be written through it.
        view = first(each_record(file; copy = false))
        @test view == copy(view) && hash(view) == hash(copy(view))
        @test HTSLib.data(view) isa ImmutableMemoryView{UInt8} && HTSLib.data(copy1) isa ImmutableMemoryView{UInt8}
        @test HTSLib.data(copy1) == HTSLib.data(view)
        @test_throws Base.CanonicalIndexError HTSLib.data(copy1)[1] = 0x00

        # `qname` is a view of the data, and `rname` the header's name, so reading either allocates nothing.
        @test HTSLib.qname(view) isa StringView{ImmutableMemoryView{UInt8}} && HTSLib.qname(view) == HTSLib.qname(copy1) == first(split(SAM_LINES[1], '\t'))
        @test HTSLib.rname(view) === HTSLib.rname(copy1) === file.header.sequence_names[1]
        @test sequence_names(file) == file.header.sequence_names && sequence_names(file) !== file.header.sequence_names
        name_allocations(record) = (@allocated HTSLib.qname(record)) + (@allocated HTSLib.rname(record))
        name_allocations(view)
        @test name_allocations(view) == 0
    end
end
