function gjk(p::Any, q::Any, dir::AbstractArray{<:Float64, 1}; tol=1e-10, max_iterations=1000)
    psimplex = support(p, dir); qsimplex = support(q, -dir); simplex = psimplex - qsimplex
    dir = -simplex
    colliding = false
    closest = [], []

    for i = 1:max_iterations
        ps = support(p, dir); qs = support(q, -dir); s = ps - qs

        if colliding
            break
        elseif s ⋅ (-dir) ≥ sum(abs2, dir)*(1.0 - tol) || any(all(simplex .== s, 1))
            λ = findcombination(simplex, -dir)
            closest = psimplex*λ, qsimplex*λ
            break
        else
            psimplex = hcat(psimplex, ps); qsimplex = hcat(qsimplex, qs); simplex = hcat(simplex, s)
            filtered, dir, colliding = findsimplex(simplex)
            psimplex = psimplex[:, filtered]; qsimplex = qsimplex[:, filtered]; simplex = simplex[:, filtered]
        end
    end
    colliding, closest
end
