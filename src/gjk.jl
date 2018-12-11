"""
    gjk(p, q, dir; atol=1e-10, max_iterations=1000)

Compute whether polytopes `p` and `q` are colliding. Provide
an initial search direction `dir` in the space of the problem.

Return result of type `Result` which informs if the computation
was successful, a collision was detected, and if not the closest
points between the two polytopes.

# Examples
```julia-repl
julia> p = [0.0 0.0 1.0; 0.0 2.0 1.0]
julia> q = [2.0 2.0 3.0; 0.0 2.0 1.0]
julia> gjk(p, q, [1.0, 0.0])
GJK.Result{Array{Float64,2}}(true, false, [1.0 2.0; 1.0 1.0])
```
"""
function gjk(p::Any, q::Any, init_dir::AbstractVector)
    ps = support(p, init_dir)
    qs = support(q, -init_dir)
    psimplex, qsimplex, dir = init_gjk(ps, qs)

    # Check if we got lucky
    if isapprox(sum(abs2, dir), 0.0)
        return true, psimplex, qsimplex
    end

    recursive_gjk(p, q, psimplex, qsimplex, dir)
end

function init_gjk(ps::SVector{N}, qs::SVector{N}) where {N}
    return SMatrix{N,1}(ps), SMatrix{N,1}(qs), qs - ps
end

function recursive_gjk(p::Any, q::Any, psimplex::SMatrix, qsimplex::SMatrix, dir::SVector)
    ps = support(p, dir)
    qs = support(q, -dir)
    if (ps - qs) ⋅ (-dir) ≥ sum(abs2, dir) || check_degeneracy(psimplex - qsimplex, ps - qs)
        pclose, qclose = closestpoints(psimplex, qsimplex, -dir)
        return false, norm(pclose - qclose), pclose, qclose
    end
    psimplex = hcat(psimplex, ps)
    qsimplex = hcat(qsimplex, qs)
    npsimplex, nqsimplex, ndir, collision = findsimplex(psimplex, qsimplex)
    if collision
        return true, 0.0, ps, qs
    else
        return recursive_gjk(p, q, npsimplex, nqsimplex, ndir)
    end
end

@generated function check_degeneracy(simplex::SMatrix{N, M, T}, s::SVector{N, T}) where {N, M, T}
    exprs = :(false)
    for i = 1:M
        releps = :(eps($T)*(1 + sum(abs2, simplex[:, $i]) + sum(abs2, s)))
        exprs = :($exprs || sum(abs2, simplex[:, $i] - s) ≤ $releps)
    end
    return exprs
end
