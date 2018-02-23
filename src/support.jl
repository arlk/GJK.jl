function support(vertices::AbstractArray{T, 2}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat}
    vertices[:, indmax(vertices'*dir)]
end

using GeometryTypes: HyperRectangle
function support(rect::HyperRectangle{N, T}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat, N}
    normvec = [if x > 0 1.0 else -1.0 end for x in normalize(dir./rect.widths, Inf)]
    rect.widths.*normvec/2.0 + rect.origin
end

using GeometryTypes: HyperCube
function support(cube::HyperCube{N, T}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat, N}
    normvec = [if x > 0 1.0 else -1.0 end for x in normalize(dir, Inf)]
    cube.width.*normvec/2.0 + cube.origin
end

using GeometryTypes: HyperSphere
function support(sphere::HyperSphere{N, T}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat, N}
    sphere.center + sphere.r*normalize(dir, 2)
end
