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
        elseif i == max_iterations
            warn("GJK has not terminated in $max_iterations iterations.")
        end
    end

    result
end
