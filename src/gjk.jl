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
    psimplex = insertcolumn(ps)
    qsimplex = insertcolumn(qs)
    dir = qs - ps

    collision, psimplex, qsimplex, dir, sz = _gjk(p, q, psimplex, qsimplex, dir)

    if collision
        return true, 0.0, ps, qs
    else
        pclose, qclose = closestpoints(psimplex, qsimplex, -dir, sz)
        return false, norm(pclose - qclose), pclose, qclose
    end
end

function _gjk(p::Any, q::Any, psimplex::SMatrix, qsimplex::SMatrix, dir::SVector)
    sz = 1
    max_iter = 100
    collision = false

    if sum(abs2, dir) ≤ 2*eps(Float64)
        return true, psimplex, qsimplex, dir, sz
    end

    ps = support(p, dir)
    qs = support(q, -dir)
    while check_gjk(psimplex - qsimplex, ps - qs, dir, sz, max_iter)
        sz += 1
        psimplex = insertcolumn(psimplex, ps, sz)
        qsimplex = insertcolumn(qsimplex, qs, sz)
        psimplex, qsimplex, dir, collision, sz = findsimplex(psimplex, qsimplex, sz)
        if collision
            break
        else
            ps = support(p, dir)
            qs = support(q, -dir)
        end

        max_iter -= 1
    end

    return collision, psimplex, qsimplex, dir, sz
end

function check_gjk(simplex::SMatrix, s::SVector, dir::SVector, sz::Int, iter::Int)
    if iter < 0
        error("GJK failed: maximum number of iterations reached")
    elseif s ⋅ (-dir) ≥ sum(abs2, dir) || check_degeneracy(simplex, s, sz)
        return false
    else
        return true
    end
end

@generated function check_degeneracy(simplex::SMatrix{N, M, T}, s::SVector{N, T}, sz::Vararg{Int, 1}) where {N, M, T}
    exprs = :(false)
    for i = 1:M
        exprs = :($i > sz[1] ? $exprs : ($exprs || sum(abs2, simplex[:, $i] - s) ≤ 2*eps($T)))
    end
    return exprs
end
