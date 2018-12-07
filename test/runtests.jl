using GJK
using Tests

@testset "GJK" begin
    include("support.jl")
    include("simplex.jl")
    include("gjk.jl")
end
