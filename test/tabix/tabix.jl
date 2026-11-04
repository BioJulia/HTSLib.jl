# Tabix Test Suite
# ================

@testset "Tabix" begin
    include("index.jl")
    include("records.jl")
    include("file.jl")
    include("iterators.jl")
end
