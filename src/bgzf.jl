# BGZF Compression
# ================

"""
    bgzip(source::AbstractString, destination::AbstractString = source * ".gz"; level = nothing, threads = 0) -> String

Compress the file `source` into the BGZF file `destination` with htslib, as `bgzip` does, and return `destination`.
Unlike `bgzip`, `source` is kept.

`level` is the compression level, from 0 to 9, or htslib's default; `threads` is the number of threads that compress blocks, or none.
"""
function bgzip(source::AbstractString, destination::AbstractString = source * ".gz"; level::Union{Nothing, Integer} = nothing, threads::Integer = 0)
    level === nothing || 0 <= level <= 9 || throw(ArgumentError("the compression level is $(level), but must be from 0 to 9"))
    file = Bindings.bgzf_open(destination, level === nothing ? "w" : "w$(level)")
    file == C_NULL && error("htslib could not open $(destination) for writing")
    try
        if threads > 0
            Bindings.bgzf_mt(file, threads, 256) == 0 || error("htslib could not start $(threads) threads to compress $(destination)")
        end
        buffer = Vector{UInt8}(undef, 1 << 20)
        open(source) do io
            while !eof(io)
                count = readbytes!(io, buffer)
                Bindings.bgzf_write(file, buffer, count) == count || error("htslib could not write to $(destination)")
            end
        end
    finally
        status = Bindings.bgzf_close(file)
        status == 0 || error("htslib could not close $(destination)")
    end
    return destination
end

"""
    bgzip_decompress(source::AbstractString, destination::AbstractString = (endswith(source, ".gz") ? source[1:(end - 3)] : source * ".decompressed"); threads = 0) -> String

Decompress the BGZF-compressed file `source` into `destination` with htslib, and return `destination`.
`threads` is the number of decompression threads, or 0.
"""
function bgzip_decompress(source::AbstractString, destination::AbstractString = (endswith(source, ".gz") ? source[1:(end - 3)] : source * ".decompressed"); threads::Integer = 0)
    file = Bindings.bgzf_open(source, "r")
    file == C_NULL && error("htslib could not open $(source) for reading")
    try
        if threads > 0
            Bindings.bgzf_mt(file, threads, 256) == 0 || error("htslib could not start $(threads) threads to decompress $(source)")
        end
        buffer = Vector{UInt8}(undef, 1 << 20)
        open(destination, "w") do io
            while true
                count = Bindings.bgzf_read(file, buffer, length(buffer))
                count < 0 && error("htslib could not read from $(source)")
                count == 0 && break
                write(io, @view buffer[1:count])
            end
        end
    finally
        status = Bindings.bgzf_close(file)
        status == 0 || error("htslib could not close $(source)")
    end
    return destination
end
