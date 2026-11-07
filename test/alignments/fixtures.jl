# Shared Alignment Test Fixtures & Utilities
# ==========================================

const SAM_HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n"

# Positions are 1-based: `read1` covers bases 100 to 119, `read2` 200 to 224 with a deletion, and `read3` 50 to 69.
const SAM_LINES = [
    "read1\t0\tchr1\t100\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII",
    "read2\t16\tchr1\t200\t30\t10M5D10M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*",
    "read3\t0\tchr2\t50\t60\t5S20M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTA\t*",
    "placed\t4\tchr2\t60\t0\t*\t*\t0\t0\tACGT\t*",
    "unplaced\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*",
]

# Convert the SAM file `source` with htslib, as `samtools view -b` does with `mode = "wb"`, and `samtools view -C -T reference` with `mode = "wc"`.
function convert_alignments(source, destination, mode; reference = nothing)
    Bindings = HTSLib.Bindings
    input = Bindings.hts_open(source, "r")
    header = Bindings.sam_hdr_read(input)
    output = Bindings.hts_open(destination, mode)
    reference === nothing || @assert Bindings.hts_set_fai_filename(output, reference) == 0
    @assert Bindings.sam_hdr_write(output, header) == 0
    record = Bindings.bam_init1()
    status = Bindings.sam_read1(input, header, record)
    while status >= 0
        @assert Bindings.sam_write1(output, header, record) >= 0
        status = Bindings.sam_read1(input, header, record)
    end
    # -1 is the end of the file, and less than -1 an invalid record.
    @assert status == -1
    Bindings.bam_destroy1(record)
    Bindings.sam_hdr_destroy(header)
    @assert Bindings.hts_close(output) == 0
    @assert Bindings.hts_close(input) == 0
    return destination
end

# Index `path` with htslib, as `samtools index` does: a BAI index, a CSI index with a positive `min_shift`, or a CRAI index of a CRAM file.
function index_alignments(path; min_shift = 0)
    @assert HTSLib.Bindings.sam_index_build3(path, C_NULL, min_shift, 0) == 0
    return path * (endswith(path, ".cram") ? ".crai" : min_shift > 0 ? ".csi" : ".bai")
end

names_of(records) = HTSLib.qname.(records)
