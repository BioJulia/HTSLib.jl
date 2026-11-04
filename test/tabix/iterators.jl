# Tabix Iterator Tests
# ====================

@testset "Tabix Iterators" begin
    directory = mktempdir()
    compressed = bgzip(write_lines(directory, "file.bed", BED_LINES; header = ["# header"]))
    tabix_index(compressed)

    TabixFile(compressed; threads = 2) do file
        # Views give the 1-based index of the sequence, the 0-based, half-open interval, and the line.
        @test [(record.sequence_index, record.start, record.stop, String(record)) for record in each_overlap(file, "chr2"; copy = false)] == [(2, 0, 5, BED_LINES[3])]
        # Reading every record skips the header, and opens the file again, so it may interleave with queries.
        overlaps_iter = each_overlap(file, "chr1"; copy = false)
        sequential = each_record(file; copy = false)
        @test String(first(iterate(overlaps_iter))) == BED_LINES[1]
        @test String(first(iterate(sequential))) == BED_LINES[1]
        @test String(first(iterate(overlaps_iter))) == BED_LINES[2]
        @test [String(record) for record in sequential] == BED_LINES[2:3]
        @test [String(record) for record in records(file)] == BED_LINES

        # Queries of the same file interleave, each reading through a handle of its own: a query that another runs within reads on from where it was, and its view stays valid until it moves on.
        chr1 = each_overlap(file, "chr1"; copy = false)
        view = first(chr1)
        @test [String(record) for record in each_overlap(file, "chr2"; copy = false)] == BED_LINES[3:3]
        @test String(view) == BED_LINES[1]
        @test [String(record) for record in chr1] == BED_LINES[2:2]

        # An unfinished query does not disturb the next.
        for record in each_overlap(file, "chr1")
            break
        end
        @test length(overlaps(file, "chr1")) == 2

        # Iterators yield records of their own, or views with `copy = false`, which `collect` copies.
        @test eltype(each_overlap(file, "chr1")) == eltype(each_record(file)) == TabixRecord{HTSLib.Owned}
        @test eltype(each_overlap(file, "chr1"; copy = false)) == eltype(each_record(file; copy = false)) == TabixRecordView
        @test all(record -> record isa TabixRecord{HTSLib.Owned}, each_overlap(file, "chr1"))
        @test collect(each_overlap(file, "chr1"; copy = false)) == overlaps(file, "chr1")
        @test collect(each_record(file; copy = false)) == records(file)
    end

    # htslib rejects a line whose sequence the index lacks, rather than yield a record without one, so every record's sequence is one of `sequence_names(file)`.
    mismatched = bgzip(write_lines(directory, "mismatched.bed", ["chr1\t10\t20", "chrX\t30\t40"]))
    # htslib also warns that the index is older than the file.
    TabixFile(mismatched; index = compressed * ".tbi") do file
        @test String(first(each_record(file))) == "chr1\t10\t20"
        @test_throws "or its sequence isn't in the index" records(file)
        @test isempty(overlaps(file, "chrX"))
    end

    # Views of lines that grow longer, so that htslib reallocates its buffer, stay right, and allocate nothing for each.
    lines = ["chr1\t$(10k)\t$(10k + 5)\t" * "x"^(4k) for k in 1:200]
    growing = bgzip(write_lines(directory, "growing.bed", lines))
    TabixFile(growing; index = tabix_index(growing)) do file
        @test [String(record) for record in each_record(file; copy = false)] == lines
        @test [String(record) for record in each_overlap(file, "chr1"; copy = false)] == lines
        total(records) = sum(length(HTSLib.line(record)) for record in records)
        total(each_record(file; copy = false))
        few = @allocated total(Iterators.take(each_record(file; copy = false), 10))
        all = @allocated total(each_record(file; copy = false))
        # The buffer is wrapped again only when htslib reallocates it, a few times as the lines grow.
        @test all - few < 100 * 32
    end
end
