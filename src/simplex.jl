function proj(u::AbstractArray{T,1}, v::AbstractArray{T,1}) where {T<:AbstractFloat}
    (u ⋅ v)/(u ⋅ u)*u
end

function findsimplex(simplex::AbstractArray{<:AbstractFloat, 2})
    if size(simplex, 2) == 2
        findline(simplex, [1, 2])
    elseif size(simplex, 2) == 3
        findtriangle(simplex, [1, 2, 3])
    elseif size(simplex, 2) == 4
        findtetrahedron(simplex, [1, 2, 3, 4])
    end
end

function findline(simplex::AbstractArray{<:AbstractFloat, 2}, idx::Array{<:Signed, 1})
    AB = simplex[:, 2] - simplex[:, 1]
    AO = simplex[:, 2]
    if AB ⋅ AO > 0
        dir = proj(AB, AO) - AO
    else
        idx = idx[[2]]
        dir = -AO
    end
    collision = all(dir .== 0.0)
    idx, dir, collision
end

function findtriangle(simplex::AbstractArray{<:AbstractFloat, 2}, idx::Array{<:Signed, 1})
    AB = simplex[:, 3] - simplex[:, 2]
    AC = simplex[:, 3] - simplex[:, 1]
    BC = simplex[:, 2] - simplex[:, 1]
    AO = simplex[:, 3]
    if (AC ⋅ AB * BC - AC ⋅ BC * AB) ⋅ AO > 0
        if AC ⋅ AO > 0
            idx = idx[[1, 3]]
            dir = proj(AC, AO) - AO
            collision = all(dir .== 0.0)
            idx, dir, collision
        else
            findline(simplex[:, [2, 3]], idx[[2, 3]])
        end
    elseif (AB ⋅ BC * AB - AB ⋅ AB * BC) ⋅ AO > 0
        findline(simplex[:, [2, 3]], idx[[2, 3]])
    else
        if size(simplex, 1) == 2
            idx, AO, true
        else
            ABC = AB × BC
            if ABC ⋅ AO > 0
                dir = -proj(ABC, AO)
            else
                idx = idx[[2, 1, 3]]
                dir = -proj(-ABC, AO)
            end
            collision = all(dir .== 0.0)
            idx, dir, collision
        end
    end
end

function findtetrahedron(simplex::AbstractArray{<:AbstractFloat, 2}, idx::Array{<:Signed, 1})
    AB = simplex[:, 4] - simplex[:, 3]
    AC = simplex[:, 4] - simplex[:, 2]
    AD = simplex[:, 4] - simplex[:, 1]
    BC = simplex[:, 3] - simplex[:, 2]
    BD = simplex[:, 3] - simplex[:, 1]
    AO = simplex[:, 4]
    BO = simplex[:, 3]
    if (AB × AC) ⋅ AO > 0
        findtriangle(simplex[:, [2, 3, 4]], idx[[2, 3, 4]])
    elseif (AC × AD) ⋅ AO > 0
        findtriangle(simplex[:, [1, 2, 4]], idx[[1, 2, 4]])
    elseif (AD × AB) ⋅ AO > 0
        findtriangle(simplex[:, [3, 1, 4]], idx[[3, 1, 4]])
    else
        idx, AO, true
    end
end

function findcombination(simplex::AbstractArray{T, 2}, vec::AbstractArray{T, 1}) where {T<:AbstractFloat}
    if size(simplex, 2) == 2
        findlinecombination(simplex, vec)
    elseif size(simplex, 2) == 3
        findtrianglecombination(simplex, vec)
    end
end

function findlinecombination(simplex::AbstractArray{T, 2}, vec::AbstractArray{T, 1}) where {T<:AbstractFloat}
    AV = simplex[:, 2] - vec
    AB = simplex[:, 2] - simplex[:, 1]
    λ = (AB ⋅ AV)/(AB ⋅ AB)
    λ = λ < 0 ? 0 : λ
    [λ; 1.0 - λ]
end

function findtrianglecombination(simplex::AbstractArray{T, 2}, vec::AbstractArray{T, 1}) where {T<:AbstractFloat}
    AO = simplex[:, 3]
    AV = simplex[:, 3] - vec
    AB = simplex[:, 3] - simplex[:, 2]
    AC = simplex[:, 3] - simplex[:, 1]
    BC = simplex[:, 2] - simplex[:, 1]

    if (AC ⋅ AB * BC - AC ⋅ BC * AB) ⋅ AV > 0
        if AC ⋅ AV > 0
            ptAC = (AO - proj(AC, AV))
            λAC = findlinecombination(simplex[:, [1, 3]], ptAC)
            [λAC[1]; 0.0; λAC[2]]
        else
            ptAB = (AO - proj(AB, AV))
            λAB = findlinecombination(simplex[:, [2, 3]], ptAB)
            [0.0; λAB[1]; λAB[2]]
        end
    elseif (AB ⋅ BC * AB - AB ⋅ AB * BC) ⋅ AV > 0
        ptAB = (AO - proj(AB, AV))
        λAB = findlinecombination(simplex[:, [2, 3]], ptAB)
        [0.0; λAB[1]; λAB[2]]
    else
        sABC = [AC AB]
        λ = (sABC' * sABC) \ (sABC' * AV)
        vcat(λ, 1-sum(λ))
    end
end
