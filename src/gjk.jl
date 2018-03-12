using StaticArrays

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
function gjk(p::Any, q::Any, dir::AbstractArray{<:Float64, 1}; atol::AbstractFloat=1e-10, max_iterations::Signed=1000)
    idx = @MMatrix zeros(3, 3)
    psimplex = @MMatrix zeros(2, 3)
    qsimplex = @MMatrix zeros(2, 3)
    simplex = @MMatrix zeros(2, 3)
    psimplex[:,1] = support(p, dir); qsimplex[:,1] = support(q, -dir);
    simplex[:] = psimplex - qsimplex
    dir[:] = -simplex[:,1]
    if isapprox(sum(abs2, dir), 0.0; atol = atol)
        return Result(true)
    end

    result = Result()
    sz = 1

    for i = 1:max_iterations
        ps = support(p, dir); qs = support(q, -dir);
        s = ps - qs

        if s ⋅ (-dir) ≥ sum(abs2, dir)*(1.0 - atol) || any(all(simplex .== s, 1))
            λ = MVector(1.,0.,0.)
            findcombination(simplex, -dir, λ, sz)
            result = Result(false, hcat(psimplex*λ, qsimplex*λ))
            break
        end

        psimplex[:,sz+1] = ps
        qsimplex[:,sz+1] = qs
        simplex[:,sz+1] = s
        idx[:] = 0
        collision = findsimplex(simplex, idx, dir, sz+1)
        sz = Int64(sum(idx))
        psimplex[:] = psimplex*idx
        qsimplex[:] = qsimplex*idx

        if collision
            result = Result(true)
            break
        elseif i ≥ max_iterations
            warn("GJK has not terminated in $max_iterations iterations.")
            break
        end
    end

    result
end
