module Collide

# package code goes here
function proj(u::AbstractArray{T,1}, v::AbstractArray{T,1}) where {T<:AbstractFloat}
    (u ⋅ v)/(u ⋅ u)*u
export gjk
end

include("support.jl")
include("simplex.jl")
include("gjk.jl")

end # module
