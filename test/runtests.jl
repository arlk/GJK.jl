using GJK
using Base.Test

@testset "GJK" begin
    include("support.jl")
    include("simplex.jl")
    include("gjk.jl")
end
