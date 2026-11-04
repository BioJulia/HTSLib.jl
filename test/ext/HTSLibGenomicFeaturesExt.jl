# GenomicFeatures Extension Tests
# ===============================
#
# With GenomicFeatures loaded, records have its interval interface, whose `group_name` is an alias of the accessor that is true to the file's specification.

import GenomicFeatures
using GenomicFeatures: group_name, has_group_name, left_position, has_left_position, right_position, has_right_position
using GenomicFeatures: GenomicInterval, GenomicIntervalCollection, strand, span, is_overlapping, STRAND_POS, STRAND_NEG, STRAND_NA

@testset "GenomicFeatures extension" begin
    directory = mktempdir()

    # Alignments: RNAME, POS and the last base of the CIGAR, all 1-based already.
    sam = write_lines(directory, "records.sam", SAM_LINES; header = [chomp(SAM_HEADER)])
    bam = convert_alignments(sam, joinpath(directory, "records.bam"), "wb")
    (read1, read2, read3, placed, unplaced) = AlignmentFile(records, bam)
    for record in (read1, read2, read3, placed)
        @test has_group_name(record) && group_name(record) == HTSLib.rname(record)
        @test has_left_position(record) && left_position(record) == HTSLib.pos(record)
        @test has_right_position(record) && right_position(record) == HTSLib.end_position(record)
    end
    @test (group_name(read1), left_position(read1), right_position(read1)) == ("chr1", 100, 119)
    @test !has_group_name(unplaced) && !has_left_position(unplaced) && !has_right_position(unplaced)
    AlignmentFile(bam) do file
        view = first(each_record(file; copy = false))
        @test (group_name(view), left_position(view), right_position(view)) == ("chr1", 100, 119)
        group_name_allocations(record) = @allocated group_name(record)
        group_name_allocations(view)
        @test group_name_allocations(view) == 0
    end

    # Tabix records: the sequence name, and the 0-based, half-open interval of htslib as 1-based and inclusive.
    compressed = bgzip(write_lines(directory, "records.bed", BED_LINES))
    TabixFile(compressed; index = tabix_index(compressed)) do file
        record = only(overlaps(file, "chr2"))
        @test has_group_name(record) && group_name(record) == HTSLib.sequence_name(record) == "chr2"
        @test (left_position(record), right_position(record)) == (1, 5)
        view = first(each_overlap(file, "chr1"; copy = false))
        @test (left_position(view), right_position(view)) == (11, 20)
        @test has_left_position(view) && has_right_position(view)
        @test has_group_name(view) && group_name(view) == "chr1"
    end

    @testset "conversions to GenomicInterval" begin
        # Tabix records: 1-based and inclusive, with an unknown strand, and the record, unchanged, as the metadata.
        bed = bgzip(write_lines(directory, "conversions.bed", ["chr1\t10\t20", "chr1\t30\t30", "chr2\t1000\t1500"]))
        TabixFile(bed; index = tabix_index(bed)) do file
            (record, empty, feature) = records(file)
            interval = GenomicInterval(record)
            @test interval isa GenomicInterval{TabixRecord{HTSLib.Owned}}
            @test (interval.groupname, interval.first, interval.last, interval.strand) == ("chr1", 11, 20, STRAND_NA)
            @test interval.metadata === record
            @test (interval.metadata.start, interval.metadata.stop) == (interval.first - 1, interval.last)
            @test GenomicInterval{TabixRecord{HTSLib.Owned}}(record).metadata === record
            @test strand(record) == STRAND_NA
            # A zero-length record, whose start is its stop in htslib, is an empty interval.
            @test (GenomicInterval(empty).first, GenomicInterval(empty).last, span(GenomicInterval(empty))) == (31, 30, 0)
            untyped = GenomicInterval{Nothing}(record)
            @test untyped isa GenomicInterval{Nothing} && untyped.metadata === nothing && (untyped.first, untyped.last) == (11, 20)
            # The interval of a view holds its owned copy, which stays valid after the next iteration.
            views = each_record(file; copy = false)
            view = first(views)
            kept = GenomicInterval(view)
            @test kept isa GenomicInterval{TabixRecord{HTSLib.Owned}} && kept.metadata == view
            @test GenomicInterval{Nothing}(view) isa GenomicInterval{Nothing}
            first(views)
            @test String(kept.metadata) == "chr1\t10\t20"
            # A collection of records, and their overlaps, are GenomicFeatures' own.
            collection = GenomicIntervalCollection{Nothing}(records(file))
            @test length(collection) == 3
            @test is_overlapping(GenomicInterval(feature), GenomicInterval("chr2", 1500, 1600))
            @test !is_overlapping(GenomicInterval(feature), GenomicInterval("chr2", 1501, 1600))
        end
        # A GFF record and a BED record of the same feature give the same interval.
        gff = bgzip(write_lines(directory, "conversions.gff", ["chr2\tsource\tgene\t1001\t1500\t.\t+\t.\tID=a"]))
        TabixFile(gff; index = tabix_index(gff)) do file
            interval = GenomicInterval{Nothing}(only(records(file)))
            @test (interval.groupname, interval.first, interval.last) == ("chr2", 1001, 1500)
        end

        # Alignments: the strand of FLAG, or unknown for an unmapped record, and none for a record without a position.
        @test (strand(read1), strand(read2), strand(placed)) == (STRAND_POS, STRAND_NEG, STRAND_NA)
        interval = GenomicInterval(read2)
        @test interval isa GenomicInterval{AlignmentRecord{HTSLib.Owned}} && interval.metadata === read2
        @test (interval.groupname, interval.first, interval.last, interval.strand) == ("chr1", 200, 224, STRAND_NEG)
        @test (GenomicInterval(placed).first, GenomicInterval(placed).last) == (60, 60)
        @test_throws "has no position" GenomicInterval(unplaced)
        @test_throws "has no position" GenomicInterval{Nothing}(unplaced)
        AlignmentFile(bam) do file
            views = each_record(file; copy = false)
            kept = GenomicInterval(first(views))
            first(views)
            @test kept isa GenomicInterval{AlignmentRecord{HTSLib.Owned}} && HTSLib.qname(kept.metadata) == "read1"
        end
    end
end
