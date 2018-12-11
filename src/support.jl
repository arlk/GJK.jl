"""
    support(cvxpoly, dir)

Compute the point of contact between a convex polytope and its
supporting hyperplane defined by the given normal direction.
"""
function support(vertices::AbstractMatrix, dir::AbstractVector)
    @inbounds vertices[:, argmax(vertices'*dir)]
end
