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
    psimplex = support(p, dir); qsimplex = support(q, -dir);
    simplex = psimplex - qsimplex
    dir = -simplex
    if isapprox(sum(abs2, dir), 0.0; atol = atol)
        return Result(true)
    end

    result = Result()

    for i = 1:max_iterations
        ps = support(p, dir); qs = support(q, -dir);
        s = ps - qs

        if s ⋅ (-dir) ≥ sum(abs2, dir)*(1.0 - atol) || any(all(simplex .== s, 1))
            λ = size(simplex, 2) == 1 || findcombination(simplex, -dir)
            result = Result(false, hcat(psimplex*λ, qsimplex*λ))
            break
        end

        psimplex = hcat(psimplex, ps); qsimplex = hcat(qsimplex, qs); simplex = hcat(simplex, s)
        filtered, dir, collision = findsimplex(simplex; atol = atol)
        psimplex = psimplex[:, filtered]; qsimplex = qsimplex[:, filtered];
        simplex = simplex[:, filtered]

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
