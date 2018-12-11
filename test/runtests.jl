using GJK
using Test
using LinearAlgebra
using StaticArrays

@testset "GJK" begin
    include("support.jl")
    #  include("simplex.jl")
    include("gjk.jl")
end
