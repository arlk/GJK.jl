"""
    proj(u, v)

Compute the vector projection of vector `v` onto vector `u`.
"""
function proj(u::AbstractVector, v::AbstractVector)
    (u ⋅ v)/(u ⋅ u)*u
end

"""
    findsimplex(simplex)

Compute the search direction for a given simplex. Return filtered
indices of the simplex to keep. Return a collision flag true if
the origin was enclosed by the simplex.
"""
function findsimplex(psimplex::SMatrix{N, 2}, qsimplex::SMatrix{N, 2}) where {N}
    simplex = psimplex - qsimplex
    AB = simplex[:, 1] - simplex[:, 2]
    AO = -simplex[:, 2]
    if AB ⋅ AO > 0
        dir = AO - proj(AB, AO)
        collision = isapprox(sum(abs2, dir), 0.0)
        return psimplex, qsimplex, dir, collision
    else
        dir = AO
        collision = isapprox(sum(abs2, dir), 0.0)
        return SMatrix{N, 1}(psimplex[:, 2]), SMatrix{N, 1}(qsimplex[:, 2]), dir, collision
    end
end

function findsimplex(psimplex::SMatrix{N, 3}, qsimplex::SMatrix{N, 3}) where {N}
    simplex = psimplex - qsimplex
    AB = simplex[:, 2] - simplex[:, 3]
    AC = simplex[:, 1] - simplex[:, 3]
    BC = simplex[:, 1] - simplex[:, 2]
    AO = -simplex[:, 3]
    if (AC ⋅ AB * BC - AC ⋅ BC * AB) ⋅ AO > 0
        if AC ⋅ AO > 0
            idx = SMatrix{3, 2}(1, 0, 0, 0, 0, 1)
            dir = AO - proj(AC, AO)
            collision = isapprox(sum(abs2, dir), 0.0)
            return psimplex*idx, qsimplex*idx, dir, collision
        else
            idx = SMatrix{3, 2}(0, 0, 1, 0, 0, 1)
            return findsimplex(psimplex*idx, qsimplex*idx)
        end
    elseif (AB ⋅ BC * AB - AB ⋅ AB * BC) ⋅ AO > 0
        idx = SMatrix{3, 2}(0, 0, 1, 0, 0, 1)
        return findsimplex(psimplex*idx, qsimplex*idx)
    else
        if isapprox(sum(abs2, AO), 0.0)
            return psimplex, qsimplex, AO, true
        elseif isapprox(sum(abs2, AC - proj(AB, AC)), 0.0)
            idx = SMatrix{3, 2}(0, 0, 1, 0, 0, 1)
            return findsimplex(psimplex*idx, qsimplex*idx)
        elseif N == 2
            return psimplex, qsimplex, AO, true
        else
            ABC = AB × BC
            dir = proj(ABC, AO)
            if ABC ⋅ AO > 0
                idx = SMatrix{3, 3}(0, 1, 0, 1, 0, 0, 0, 0, 1)
                collision = isapprox(sum(abs2, dir), 0.0)
                return psimplex*idx, qsimplex*idx, dir, collision
            else
                collision = isapprox(sum(abs2, dir), 0.0)
                return psimplex, qsimplex, dir, collision
            end
        end
    end
end

function findsimplex(psimplex::SMatrix{N, 4}, qsimplex::SMatrix{N, 4}) where {N}
    simplex = psimplex - qsimplex
    AB = simplex[:, 3] - simplex[:, 4]
    AC = simplex[:, 2] - simplex[:, 4]
    AD = simplex[:, 1] - simplex[:, 4]
    BC = simplex[:, 2] - simplex[:, 3]
    BD = simplex[:, 1] - simplex[:, 3]
    AO = -simplex[:, 4]
    if (AB × AC) ⋅ AO < 0
        idx = SMatrix{4, 3}(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
        return findsimplex(psimplex*idx, qsimplex*idx)
    elseif (AD × AB) ⋅ AO < 0
        idx = SMatrix{4, 3}(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        return findsimplex(psimplex*idx, qsimplex*idx)
    elseif (AC × AD) ⋅ AO < 0
        idx = @SMatrix zeros(M, M-1)
        idx = SMatrix{4, 3}(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1)
        return findsimplex(psimplex*idx, qsimplex*idx)
    else
        return psimplex, qsimplex, AO, true
    end
end

"""
    findcombination(psimplex, qsimplex, vec)

Compute the coefficients of a convex combination of the given
simplex vertices that equals the vector `vec`. Project `vec`
onto the closest vertex or edge of the simplex if no solution
is available.
"""
function closestpoints(psimplex::SMatrix{N, M}, qsimplex::SMatrix{N, M}, dir2origin::SVector{N}) where {N, M}
    λ = SVector{M}(cvxcombination(psimplex - qsimplex, dir2origin))
    return psimplex*λ, qsimplex*λ
end

"""
    findlinecombination(simplex, vec)

Called by `findcombination()` to solve the convex combination
problem for a line simplex defined by two points.
"""
function cvxcombination(simplex::SMatrix{N, 1}, vec::SVector{N}) where {N}
    return 1.0
end

function cvxcombination(simplex::SMatrix{N, 2}, vec::SVector{N}) where {N}
    AV = vec - simplex[:, 2]
    AB = simplex[:, 1] - simplex[:, 2]
    λ = (AB ⋅ AV)/(AB ⋅ AB)
    if λ < 0
        return 0.0, 1.0 - λ
    else
        return λ, 1.0 - λ
    end
end

function cvxcombination(simplex::SMatrix{N, 3}, vec::SVector{N}) where {N}
    AO = -simplex[:, 3]
    AV = vec - simplex[:, 3]
    AB = simplex[:, 2] - simplex[:, 3]
    AC = simplex[:, 1] - simplex[:, 3]
    BC = simplex[:, 1] - simplex[:, 2]

    if (AC ⋅ AB * BC - AC ⋅ BC * AB) ⋅ AV > 0
        if AC ⋅ AV > 0
            dir = AO - proj(AC, AV)
            idx = SMatrix{3, 2}(1, 0, 0, 0, 0, 1)
            λAC1, λAC2 = cvxcombination(simplex*idx, -dir)
            return λAC1, 0.0, λAC2
        else
            dir = AO - proj(AB, AV)
            idx = SMatrix{3, 2}(0, 0, 1, 0, 0, 1)
            λAB1, λAB2 = cvxcombination(simplex*idx, -dir)
            return 0.0, λAB1, λAB2
        end
    elseif (AB ⋅ BC * AB - AB ⋅ AB * BC) ⋅ AV > 0
        dir = AO - proj(AB, AV)
        idx = SMatrix{3, 2}(0, 0, 1, 0, 0, 1)
        λAB1, λAB2 = cvxcombination(simplex*idx, -dir)
        return 0.0, λAB1, λAB2
    else
        sABC = [AC AB]
        λ = (sABC' * sABC) \ (sABC' * AV)
        return λ[1], λ[2], 1-sum(λ)
    end
end
