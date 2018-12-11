module GJK

export gjk

using LinearAlgebra
using StaticArrays

struct Result{T<:AbstractArray{<:AbstractFloat, 2}}
    status::Bool
    collision::Bool
    contact::T
end

Result(dir::AbstractVector) = Result(false, false, Array{Float64, 2}(undef, length(dir), 2))
Result(res::Bool, dir::AbstractVector) = Result(true, res, Array{Float64, 2}(undef, length(dir), 2))
Result(res::Bool, arr::AbstractArray{<:AbstractFloat, 2}) = Result(true, res, arr)

include("support.jl")
include("simplex.jl")
include("gjk.jl")

end # module
