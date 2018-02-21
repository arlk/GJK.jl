module GJK

export gjk

struct Result{T<:AbstractArray{<:AbstractFloat, 2}}
    status::Bool
    collision::Bool
    contact::T
end

Result() = Result(false, false, Array{Float64, 2}(0, 0))
Result(res::Bool) = Result(true, res, Array{Float64, 2}(0, 0))
Result(res::Bool, arr::AbstractArray{<:AbstractFloat, 2}) = Result(true, res, arr)

include("support.jl")
include("simplex.jl")
include("gjk.jl")

end # module
