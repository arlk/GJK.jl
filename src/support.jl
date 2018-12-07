"""
    support(cvxpoly, dir)

Compute the point of contact between a convex polytope and its
supporting hyperplane defined by the given normal direction.
"""
function support(vertices::AbstractArray{T, 2}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat}
    vertices[:, argmax(vertices'*dir)]
end
