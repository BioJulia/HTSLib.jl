# What HTSLib.Bindings binds of htslib, and how its names and types differ from what Clang.jl generates.
# `generate.jl` includes this file; see README.md.

# Rules
# -----
#
# Clang.jl's argument and return types are rewritten by these rules, then by each entry's overrides.
# Arguments and returns whose types are enums become `Cint`.

# Strings are `Cstring`, which converts a `String` to a pointer to its NUL-terminated bytes.
# Optional strings, which take `C_NULL`, are overridden as `Ptr{UInt8}`.
const ARGUMENT_TYPES = Dict(
    :(Ptr{Cchar}) => :Cstring,
    :(Ptr{kstring_t}) => :(Ref{kstring_t}),
    :(Ptr{Cint}) => :(Ref{Cint}),
    :(Ptr{hts_pos_t}) => :(Ref{hts_pos_t}),
)

const RETURN_TYPES = Dict(
    :(Ptr{Cchar}) => :Cstring,
)

# Argument names that would clash with Julia's keywords.
const ARGUMENT_NAMES = Dict(
    :_end => :end_,
)

# Functions
# ---------
#
# Each entry is a C function, in output order under its section, alone or with overrides:
# - `arguments`, the types of arguments by their C names;
# - `returns`, the return type;
# - `wrap`, a function applied to the result;
# - `note`, a comment after the C declaration.

const FUNCTIONS = [
    "Files" => [
        :hts_version => (wrap = :unsafe_string,),
        :hts_set_log_level,
        :hts_get_log_level,
        :hts_open,
        :hts_close,
        :hts_get_format,
        :hts_set_fai_filename,
        :hts_idx_destroy,
        :hts_set_threads,
        :hts_get_bgzfp,
    ],
    "BGZF" => [
        :bgzf_open,
        :bgzf_close,
        :bgzf_mt,
        :bgzf_write => (arguments = (data = :(Ptr{UInt8}),),),
        :bgzf_read => (arguments = (data = :(Ptr{UInt8}),),),
        :bgzf_peek,
        :bgzf_getline,
    ],
    "Tabix" => [
        :tbx_index_build3 => (arguments = (fnidx = :(Ptr{UInt8}), conf = :(Ref{tbx_conf_t})),),
        :tbx_index_load3 => (arguments = (fnidx = :(Ptr{UInt8}),),),
        :tbx_destroy,
        :tbx_seqnames => (returns = :(Ptr{Cstring}), note = "free the array but not the values"),
        :tbx_name2id,
        :tbx_readrec => (arguments = (tbxv = :(Ptr{tbx_t}), sv = :(Ref{kstring_t})),),
        :hts_itr_destroy,
    ],
    "SAM, BAM and CRAM" => [
        :sam_hdr_read,
        :sam_hdr_write,
        :sam_hdr_destroy,
        :sam_hdr_dup,
        :sam_hdr_parse => (arguments = (text = :(Ptr{UInt8}),),),
        :sam_hdr_length,
        :sam_hdr_str => (returns = :(Ptr{UInt8}),),
        :sam_hdr_nref,
        :sam_hdr_tid2name,
        :sam_hdr_tid2len,
        :sam_hdr_name2tid,
        :bam_init1,
        :bam_destroy1,
        :bam_dup1,
        :bam_endpos,
        :sam_read1,
        :sam_write1,
        :sam_parse1,
        :sam_format1,
        :sam_index_build3 => (arguments = (fnidx = :(Ptr{UInt8}),),),
        :sam_index_load3 => (arguments = (fnidx = :(Ptr{UInt8}),),),
        :sam_itr_queryi,
        :sam_itr_querys,
    ],
    "Auxiliary fields" => [
        :bam_aux_get => (arguments = (tag = :(Ptr{UInt8}),),),
        :bam_aux2i,
        :bam_aux2f,
        :bam_aux2A,
        :bam_aux2Z,
        :bam_auxB_len,
        :bam_auxB2i,
        :bam_auxB2f,
    ],
]

# Structs
# -------
#
# Structs without bit fields or unions, which Clang.jl writes field for field.
# Each entry has a docstring, and may override the types of fields by their C names with `fields`.

const STRUCTS = [
    :kstring_t => (
        doc = "A growable string from htslib/kstring.h: `l` bytes are used of the `m` allocated at `s`, which htslib allocates and the caller frees with `Libc.free`.",
        fields = (s = :(Ptr{UInt8}),),
    ),
    :bam1_core_t => (
        doc = """
        The fixed-length fields of an alignment, from htslib/sam.h.
        Positions are 0-based, and `tid` and `mtid` are 0-based indexes of reference sequences, or -1.""",
    ),
    :tbx_conf_t => (
        doc = """
        The columns of a tabix-indexed file, from htslib/tbx.h.
        Columns are 1-based, and `preset` is one of `TBX_GENERIC`, `TBX_SAM` and `TBX_VCF`, with `TBX_UCSC` set for 0-based, half-open positions.""",
    ),
    :tbx_t => (doc = "A tabix index, from htslib/tbx.h.",),
]

# Constants
# ---------
#
# Macros and enum members, by name and the name of the type of their values, under a comment.

const CONSTANTS = [
    "Presets of tbx_conf_t." => [
        :TBX_GENERIC => :Int32,
        :TBX_SAM => :Int32,
        :TBX_VCF => :Int32,
        :TBX_UCSC => :Int32,
    ],
    "Values of enum htsExactFormat." => [
        :sam => :Cint,
        :bam => :Cint,
        :cram => :Cint,
        :vcf => :Cint,
        :bcf => :Cint,
    ],
    "Values of enum hts_fmt_option." => [
        :CRAM_OPT_DECODE_MD => :Cint,
    ],
    "Flags of tbx_index_load3 and sam_index_load3." => [
        :HTS_IDX_SAVE_REMOTE => :Cint,
        :HTS_IDX_SILENT_FAIL => :Cint,
    ],
    "enum htsLogLevel, from htslib/hts_log.h." => [
        :HTS_LOG_OFF => :Cint,
        :HTS_LOG_ERROR => :Cint,
        :HTS_LOG_WARNING => :Cint,
        :HTS_LOG_INFO => :Cint,
        :HTS_LOG_DEBUG => :Cint,
        :HTS_LOG_TRACE => :Cint,
    ],
]

# Layouts
# -------
#
# The structs that Bindings writes by hand, because they have bit fields or unions, which are checked against htslib.
# Each entry maps the fields of Bindings' structs to the C fields at the same offsets, where their names differ, with paths such as `"version.major"` into nested structs.
# A `prefix` struct has only the leading fields, so its size isn't checked.

const LAYOUTS = [
    :htsFile => (prefix = true, fields = (flags = :is_bin,)),
    :htsFormat => (prefix = false, fields = (version_major = "version.major", version_minor = "version.minor")),
    :bam1_t => (prefix = false, fields = (;)),
    :hts_itr_t => (prefix = true, fields = (flags = :read_rest,)),
]
