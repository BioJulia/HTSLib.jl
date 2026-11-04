# Alignment Test Suite
# ====================

@testset "Alignments" begin
    include("fixtures.jl")
    include("records.jl")
    include("formats/sam.jl")
    include("formats/bam.jl")
    include("formats/cram.jl")
    include("iterators.jl")
end
