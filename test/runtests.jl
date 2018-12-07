using GJK
using Test
using LinearAlgebra

@testset "GJK" begin
    include("support.jl")
    include("simplex.jl")
    include("gjk.jl")
end
