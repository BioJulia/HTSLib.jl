"""
    HTSLib.Bindings

Low-level bindings to libhts from htslib_jll.

Each function calls the C function of the same name, whose declaration is in the comment above it, and takes and returns raw pointers; the caller checks errors and manages memory.
The structs mirror their C counterparts in htslib 1.10 or later, where `hts_pos_t` is a 64-bit integer, and keep the C names of their fields.
Structs that are only used through pointers are opaque, aliases of `Cvoid` named after the C type.
They cover the parts of the API that HTSLib uses, and others that packages need, such as writing alignments and building their indexes.
"""
module Bindings

# This module is assembled in six parts, from code written here by hand and from two files that gen/generate.jl writes with Clang.jl from htslib's headers:
# 1. here: `hts_pos_t` and the opaque structs, which the generated files refer to;
# 2. generated, `Bindings/types.jl`: structs without bit fields or unions, and integer constants;
# 3. here: structs with bit fields or unions, which Clang.jl can't write as Julia structs; they use structs of part 2, and the generator checks their layouts;
# 4. here: string constants and extern globals;
# 5. generated, `Bindings/functions.jl`: functions;
# 6. here: macros, and inline and variadic functions, which Clang.jl can't bind; they call functions of part 5.
#
# The parts are in this order because each uses names from the parts before it.
# gen/selection.jl lists what the generated files contain, and how their names and types differ from htslib's; gen/README.md describes how to regenerate them.
# Don't edit the generated files: edit gen/selection.jl and run the generator.

import htslib_jll

# htslib_jll defines `libhts` only on the platforms it builds htslib for, which exclude Windows.
# Elsewhere the bindings still load, so that the package does, and each call fails; `HTSLib.is_available()` says which.
const libhts = htslib_jll.is_available() ? htslib_jll.libhts : "libhts"

# Part 1: types that the generated files refer to
# -----------------------------------------------

const hts_pos_t = Int64

# Opaque structs, used only through pointers.
const BGZF = Cvoid
const hts_idx_t = Cvoid
const sam_hdr_t = Cvoid

# Part 2, generated: structs without bit fields or unions, and integer constants
# ------------------------------------------------------------------------------

include("Bindings/types.jl")

kstring_t() = kstring_t(0, 0, C_NULL)

# Part 3: structs with bit fields or unions
# -----------------------------------------
#
# Each has the fields up to the last that is read, with bit fields as their containing integer and unions as a pointer.
# gen/generate.jl checks their offsets and sizes against htslib's, as LAYOUTS in gen/selection.jl lists them.

"""
    htsFile

The leading fields of an open file, from htslib/hts.h, up to the union of its stream pointers.
Only pointers to it are used, and only these fields are read.
`samFile` is the same type.
"""
struct htsFile
    flags::UInt32 # The bit fields is_bin, is_write, is_be, is_cram, is_bgzf and dummy, from the least significant bit.
    lineno::Int64
    line::kstring_t
    fn::Ptr{UInt8}
    fn_aux::Ptr{UInt8}
    fp::Ptr{Cvoid} # The union of BGZF *bgzf, struct cram_fd *cram and struct hFILE *hfile.
end

const samFile = htsFile

"""
    htsFormat

The format of a file, from htslib/hts.h.
`category`, `format` and `compression` are values of the enums `htsFormatCategory`, `htsExactFormat` and `htsCompression`.
"""
struct htsFormat
    category::Cint
    format::Cint
    version_major::Cshort
    version_minor::Cshort
    compression::Cint
    compression_level::Cshort
    specific::Ptr{Cvoid}
end

"""
    bam1_t

An alignment, from htslib/sam.h: its fixed-length fields, and `l_data` bytes at `data` that hold its name, CIGAR, sequence, qualities and auxiliary fields, as the `bam_get_*` macros locate them.
"""
struct bam1_t
    core::bam1_core_t
    id::UInt64
    data::Ptr{UInt8}
    l_data::Cint
    m_data::UInt32
    mempolicy::UInt32 # The bit field mempolicy and reserved bits.
end

"""
    hts_itr_t

The leading fields of an iterator over a region, from htslib/hts.h, up to the interval of the record last read.
Only pointers to it are used, and only these fields are read.
"""
struct hts_itr_t
    flags::UInt32 # The bit fields read_rest, finished, is_cram, nocoor, multi and dummy.
    tid::Cint
    n_off::Cint
    i::Cint
    n_reg::Cint
    beg::hts_pos_t
    var"end"::hts_pos_t
    reg_list::Ptr{Cvoid}
    curr_tid::Cint
    curr_reg::Cint
    curr_intv::Cint
    curr_beg::hts_pos_t
    curr_end::hts_pos_t
end

# Part 4: string constants and extern globals
# -------------------------------------------

# #define BAM_CIGAR_STR "MIDNSHP=XB"
const BAM_CIGAR_STR = "MIDNSHP=XB"

# extern const char seq_nt16_str[];
const seq_nt16_str = "=ACMGRSVTWYHKDBN"

# extern const tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf;
tbx_conf_gff() = unsafe_load(cglobal((:tbx_conf_gff, libhts), tbx_conf_t))
tbx_conf_bed() = unsafe_load(cglobal((:tbx_conf_bed, libhts), tbx_conf_t))
tbx_conf_psltbl() = unsafe_load(cglobal((:tbx_conf_psltbl, libhts), tbx_conf_t))
tbx_conf_sam() = unsafe_load(cglobal((:tbx_conf_sam, libhts), tbx_conf_t))
tbx_conf_vcf() = unsafe_load(cglobal((:tbx_conf_vcf, libhts), tbx_conf_t))

# Part 5, generated: functions
# ----------------------------

include("Bindings/functions.jl")

# Part 6: macros, and inline and variadic functions
# -------------------------------------------------

# int hts_set_opt(htsFile *fp, enum hts_fmt_option opt, ...);
# Bound for options whose value is an int, passed as a variadic argument.
hts_set_opt(fp, opt, value) = @ccall libhts.hts_set_opt(fp::Ptr{htsFile}, opt::Cint; value::Cint)::Cint

# `tbx_itr_queryi`, `tbx_itr_querys` and `tbx_itr_next` are macros in htslib/tbx.h, which these functions expand.

# #define tbx_itr_queryi(tbx, tid, beg, end) hts_itr_query((tbx)->idx, (tid), (beg), (end), tbx_readrec)
function tbx_itr_queryi(tbx, tid, beg, end_)
    idx = unsafe_load(tbx).idx
    return ccall((:hts_itr_query, libhts), Ptr{hts_itr_t}, (Ptr{hts_idx_t}, Cint, hts_pos_t, hts_pos_t, Ptr{Cvoid}), idx, tid, beg, end_, cglobal((:tbx_readrec, libhts)))
end

# #define tbx_itr_querys(tbx, s) hts_itr_querys((tbx)->idx, (s), (hts_name2id_f)(tbx_name2id), (tbx), hts_itr_query, tbx_readrec)
function tbx_itr_querys(tbx, s)
    idx = unsafe_load(tbx).idx
    return ccall(
        (:hts_itr_querys, libhts), Ptr{hts_itr_t},
        (Ptr{hts_idx_t}, Cstring, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        idx, s, cglobal((:tbx_name2id, libhts)), tbx, cglobal((:hts_itr_query, libhts)), cglobal((:tbx_readrec, libhts)),
    )
end

# #define tbx_itr_next(htsfp, tbx, itr, r) hts_itr_next(hts_get_bgzfp(htsfp), (itr), (r), (tbx))
function tbx_itr_next(htsfp, tbx, itr, r)
    return ccall((:hts_itr_next, libhts), Cint, (Ptr{BGZF}, Ptr{hts_itr_t}, Ref{kstring_t}, Ptr{tbx_t}), hts_get_bgzfp(htsfp), itr, r, tbx)
end

# `sam_itr_next` is a static inline function in htslib/sam.h, which this function reproduces for iterators of a single region.
# static inline int sam_itr_next(htsFile *htsfp, hts_itr_t *itr, bam1_t *r)
function sam_itr_next(htsfp, itr, r)
    file = unsafe_load(htsfp)
    is_cram = (file.flags >> 3) & 0x01 == 0x01
    is_bgzf = (file.flags >> 4) & 0x01 == 0x01
    (is_bgzf || is_cram) || return Cint(-2)
    itr == C_NULL && return Cint(-2)
    return ccall((:hts_itr_next, libhts), Cint, (Ptr{BGZF}, Ptr{hts_itr_t}, Ptr{bam1_t}, Ptr{htsFile}), is_bgzf ? file.fp : C_NULL, itr, r, htsfp)
end

end # module Bindings
