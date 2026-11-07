# Tabix Record and View Tests
# ============================

@testset "Tabix Records & Views" begin
    directory = mktempdir()
    plain = write_lines(directory, "records.bed", BED_LINES; header = ["# header"])
    compressed = bgzip(plain)
    tbi = tabix_index(compressed)

    TabixFile(compressed; index = tbi) do file
        recs = overlaps(file, "chr2")
        @test length(recs) == 1
        rec = first(recs)
        @test rec isa TabixRecord{HTSLib.Owned}
        @test HTSLib.sequence_name(rec) == "chr2"
        @test rec.start == 0
        @test rec.stop == 5
        @test String(rec) == BED_LINES[3]

        # The sequence name, after the tabix specification's "Column for the sequence name", and the positions.
        @test HTSLib.sequence_name(rec) == "chr2"
        @test HTSLib.position(rec) == 1
        @test HTSLib.interval_0based(rec) == (0, 5)

        # TabixRecordView
        view = first(each_overlap(file, "chr2"; copy = false))
        @test view isa TabixRecordView
        @test view.sequence_index == 2
        @test view.start == 0
        @test view.stop == 5
        @test String(view) == BED_LINES[3]
        @test HTSLib.position(view) == 1
        @test HTSLib.interval_0based(view) == (0, 5)
        @test HTSLib.sequence_name(view) == "chr2"

        # `copy` gives a TabixRecord, which owns its data, and is immutable, so a copy of it is itself.
        rec_copy = copy(view)
        @test rec_copy == TabixRecord("chr2", 0, 5, BED_LINES[3])
        @test copy(rec_copy) === rec_copy
        # A view equals its copy, with the same hash.
        @test view == rec_copy && hash(view) == hash(rec_copy)
        @test repr(rec_copy) == "TabixRecord(\"chr2\", 0, 5, \"chr2\\t0\\t5\\tc\")"
        @test startswith(repr(view), "TabixRecordView(\"chr2\"")

        # `line` gives the bytes of the line, which can't be written through it.
        @test HTSLib.line(view) isa ImmutableMemoryView{UInt8} && HTSLib.line(rec_copy) isa ImmutableMemoryView{UInt8}
        @test HTSLib.line(rec_copy) == codeunits(BED_LINES[3])
        @test_throws Base.CanonicalIndexError HTSLib.line(rec_copy)[1] = 0x00
    end
end
