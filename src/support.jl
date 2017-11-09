function support(vertices::AbstractArray{T, 2}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat}
    vertices[:, indmax(vertices'*dir)]
end

using GeometryTypes: HyperRectangle
function support(rect::HyperRectangle{N, T}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat, N}
    rect.widths.*normalize(dir./rect.widths, Inf)/2.0 + rect.origin
end

using GeometryTypes: HyperCube
function support(cube::HyperCube{N, T}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat, N}
    cube.width*normalize(dir, Inf)/2.0 + cube.origin
end

using GeometryTypes: HyperSphere
function support(sphere::HyperSphere{N, T}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat, N}
    sphere.center + sphere.r*dir/norm(dir)
end
