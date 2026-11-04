# Messages and errors

## htslib's messages

htslib reports what it notices on the standard error, each message with its level and the function that wrote it:

```
[W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
[E::hts_parse_region] Unexpected string "abc" after region
[E::idx_find_and_load] Could not retrieve index file for 'reads.bam'
```

`W` is a warning, about something htslib corrected or doubts, and `E` an error, about something it couldn't do.
HTSLib.jl leaves these as htslib writes them, as `tabix` and samtools do: it doesn't change htslib's log level, or ask htslib to fail without a message.
They are htslib's account of a file, which is what a reference implementation is for.

htslib writes them from C, so Julia's logging doesn't see them, though `redirect_stderr` does, since it redirects the process's standard error.
htslib's own setting, its log level, is in the low-level bindings, `HTSLib.Bindings.hts_set_log_level`, with the levels `HTS_LOG_OFF` to `HTS_LOG_TRACE`, as the `--verbosity` option of htslib's tools sets it.

Some common messages:

| Message | Meaning |
|:--|:--|
| `Coordinate <= 0 detected. Did you forget to use the -0 option?` | a start before the first base, read as the first base; often a 0-based file indexed as 1-based |
| `Failed to parse TBX_…, was wrong -p [type] used?` | a line whose columns the index can't read |
| `Unexpected string … after region` | a region string that htslib can't parse, which has no records |
| `Could not retrieve index file for …` | no index where htslib looked for one |
| `The index file is older than the data file` | an index that may not describe the file, which htslib uses anyway |
| `Region … cannot be stored in a tbi index. Try using a csi index` | a position past 2^29, which needs `csi = true` |

## What the package throws

Where htslib fails, the package throws an `ErrorException` that says what htslib couldn't do, after htslib's own message, if it wrote one:

```@example messages
using HTSLib

missing_file = joinpath(mktempdir(), "missing.bed.gz")
try
    TabixFile(missing_file)
catch error
    replace(sprint(showerror, error), dirname(missing_file) => "…")
end
```

A call that can't work throws an `ArgumentError` before calling htslib, for example:
- a query of a closed file, or of a file that isn't BGZF-compressed;
- a `tabix_index` without a preset, of a file whose preset can't be detected, or of a BAM, CRAM or BCF file;
- a preset that htslib doesn't have, or a `min_shift` without `csi = true`.

As in htslib, some things aren't errors: a query of a sequence that the file doesn't have, or of a region that htslib can't parse, has no records.
