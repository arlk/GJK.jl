"""
    proj(u, v)

Compute the vector projection of vector `v` onto vector `u`.
"""
function proj(u::AbstractArray{T,1}, v::AbstractArray{T,1}) where {T<:AbstractFloat}
    (u ⋅ v)/(u ⋅ u)*u
end

"""
    findsimplex(simplex)

Compute the search direction for a given simplex. Return filtered
indices of the simplex to keep. Return a collision flag true if
the origin was enclosed by the simplex.
"""
function findsimplex(simplex::AbstractArray{T, 2}, idx::AbstractArray{T, 2}, dir::AbstractArray{T, 1}, sz::Signed) where {T <: AbstractFloat}
    if sz == 2
        return findline(simplex, idx, dir)
    elseif sz == 3
        return findtriangle(simplex, idx, dir)
    elseif sz == 4
        return findtetrahedron(simplex, idx, dir)
    end
end

"""
    findline(simplex)

Called by `findsimplex()` to perform simplex operations
for a line simplex defined by two points.
"""
function findline(simplex::AbstractArray{T, 2}, idx::AbstractArray{T, 2}, dir::AbstractArray{T, 1}) where {T <: AbstractFloat}
    AB = simplex[:, 1] - simplex[:, 2]
    AO = -simplex[:, 2]
    if AB ⋅ AO > 0
        idx[1, 1] = 1
        idx[2, 2] = 1
        dir[:] = AO - proj(AB, AO)
    else
        idx[2, 1] = 1
        dir[:] = AO
    end
    simplex[:] = simplex*idx
    return isapprox(sum(abs2, dir), 0.0)
end

"""
    findtriangle(simplex)

Called by `findsimplex()` to perform simplex operations
for a triangle simplex defined by three points.
"""
function findtriangle(simplex::AbstractArray{T, 2}, idx::AbstractArray{T, 2}, dir::AbstractArray{T, 1}) where {T <: AbstractFloat}
    AB = simplex[:, 2] - simplex[:, 3]
    AC = simplex[:, 1] - simplex[:, 3]
    BC = simplex[:, 1] - simplex[:, 2]
    AO = -simplex[:, 3]
    if (AC ⋅ AB * BC - AC ⋅ BC * AB) ⋅ AO > 0
        if AC ⋅ AO > 0
            idx[1, 1] = 1
            idx[3, 2] = 1
            dir[:] = AO - proj(AC, AO)
            simplex[:] = simplex*idx
            return isapprox(sum(abs2, dir), 0.0)
        else
            idx[2, 1] = 1
            idx[3, 2] = 1
            simplex[:] = simplex*idx
            return findline(simplex, idx, dir)
        end
    elseif (AB ⋅ BC * AB - AB ⋅ AB * BC) ⋅ AO > 0
        idx[2, 1] = 1
        idx[3, 2] = 1
        simplex[:] = simplex*idx
        return findline(simplex, idx, dir)
    else
        if isapprox(sum(abs2, AO), 0.0)
            dir[:] = AO
            return true
        elseif isapprox(sum(abs2, AC - proj(AB, AC)), 0.0)
            idx[2, 1] = 1
            idx[3, 2] = 1
            simplex[:] = simplex*idx
            return findline(simplex, idx, dir)
        elseif size(simplex, 1) == 2
            dir[:] = AO
            return true
        else
            ABC = AB × BC
            dir[:] = proj(ABC, AO)
            if ABC ⋅ AO > 0
                idx[2, 1] = 1
                idx[1, 2] = 1
                idx[3, 3] = 1
            else
                idx[1, 1] = 1
                idx[2, 2] = 1
                idx[3, 3] = 1
            end
            return isapprox(sum(abs2, dir), 0.0)
        end
    end
end

"""
    findtetrahedron(simplex)

Called by `findsimplex()` to perform simplex operations
for a tetrahedron simplex defined by four points.
"""
function findtetrahedron(simplex::AbstractArray{T, 2}, idx::AbstractArray{T, 2}, dir::AbstractArray{T, 1}) where {T <: AbstractFloat}
    AB = simplex[:, 3] - simplex[:, 4]
    AC = simplex[:, 2] - simplex[:, 4]
    AD = simplex[:, 1] - simplex[:, 4]
    BC = simplex[:, 2] - simplex[:, 3]
    BD = simplex[:, 1] - simplex[:, 3]
    AO = -simplex[:, 4]
    if (AB × AC) ⋅ AO < 0
        idx[2, 1] = 1
        idx[3, 2] = 1
        idx[4, 3] = 1
        simplex[:] = simplex*idx
        return findtriangle(simplex, idx, dir)
    elseif (AD × AB) ⋅ AO < 0
        idx[3, 1] = 1
        idx[1, 2] = 1
        idx[4, 3] = 1
        simplex[:] = simplex*idx
        return findtriangle(simplex, idx, dir)
    elseif (AC × AD) ⋅ AO < 0
        idx[1, 1] = 1
        idx[2, 2] = 1
        idx[4, 3] = 1
        simplex[:] = simplex*idx
        return findtriangle(simplex, idx, dir)
    else
        dir[:] = AO
        return true
    end
end

"""
    findcombination(simplex, vec)

Compute the coefficients of a convex combination of the given
simplex vertices that equals the vector `vec`. Project `vec`
onto the closest vertex or edge of the simplex if no solution
is available.
"""
function findcombination(simplex::AbstractArray{T, 2}, vec::AbstractArray{T, 1}, comb::AbstractArray{T, 1}, sz::Signed) where {T<:AbstractFloat}
    if sz == 2
        findlinecombination(simplex, vec, comb)
    elseif sz == 3
        findtrianglecombination(simplex, vec, comb)
    end
end

"""
    findlinecombination(simplex, vec)

Called by `findcombination()` to solve the convex combination
problem for a line simplex defined by two points.
"""
function findlinecombination(simplex::AbstractArray{T, 2}, vec::AbstractArray{T, 1}, comb::AbstractArray{T, 1}) where {T<:AbstractFloat}
    AV = vec - simplex[:, 2]
    AB = simplex[:, 1] - simplex[:, 2]
    λ = (AB ⋅ AV)/(AB ⋅ AB)
    comb[1] = λ < 0 ? 0 : λ
    comb[2] = 1 - comb[1]
end

"""
    findtrianglecombination(simplex, vec)

Called by `findcombination()` to solve the convex combination
problem for a triangle simplex defined by three points.
"""
function findtrianglecombination(simplex::AbstractArray{T, 2}, vec::AbstractArray{T, 1}, comb::AbstractArray{T, 1}) where {T<:AbstractFloat}
    AO = -simplex[:, 3]
    AV = vec - simplex[:, 3]
    AB = simplex[:, 2] - simplex[:, 3]
    AC = simplex[:, 1] - simplex[:, 3]
    BC = simplex[:, 1] - simplex[:, 2]

    if (AC ⋅ AB * BC - AC ⋅ BC * AB) ⋅ AV > 0
        if AC ⋅ AV > 0
            dir = AO - proj(AC, AV)
            findlinecombination(simplex*SMatrix{3,2}(1,0,0,0,0,1), -dir, comb)
            comb[3] = comb[2]
            comb[2] = 0.0
        else
            dir = AO - proj(AB, AV)
            findlinecombination(simplex*SMatrix{3,2}(0,1,0,0,0,1), -dir, comb)
            comb[3] = comb[2]
            comb[2] = comb[1]
            comb[1] = 0.0
        end
    elseif (AB ⋅ BC * AB - AB ⋅ AB * BC) ⋅ AV > 0
        dir = AO - proj(AB, AV)
        findlinecombination(simplex*SMatrix{3,2}(0,1,0,0,0,1), -dir, comb)
        comb[3] = comb[2]
        comb[2] = comb[1]
        comb[1] = 0.0
    else
        sABC = [AC AB]
        λ = (sABC' * sABC) \ (sABC' * AV)
        comb[1] = λ[1]
        comb[2] = λ[2]
        comb[3] = 1-sum(λ)
    end
end
