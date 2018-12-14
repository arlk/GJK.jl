using Makie
using QHull: chull
using GJK: closest_points
using StaticArrays

function randomcvx(sz::Int, center::Array)
    pts = rand(sz, 2) .+ center
    ch = chull(pts)
    sz = length(ch.vertices)
    return randomcvx(ch.points[ch.vertices, :], Val(sz))
end

function randomcvx(hull::Array, ::Val{N}) where {N}
    return SMatrix{N, 2}(hull)
end

function drawjulialogo(scene)
    centerA = [ 0. √2.]
    centerB = [ 1.  0.]
    centerC = [-1.  0.]
    polyA = randomcvx(15, centerA)
    polyB = randomcvx(15, centerB)
    polyC = randomcvx(15, centerC)
    AB = closest_points(polyA', polyB', SVector{2}(centerA - centerB))
    BC = closest_points(polyB', polyC', SVector{2}(centerB - centerC))
    CA = closest_points(polyC', polyA', SVector{2}(centerC - centerA))
    vAB = [Point2f0(x) for x in AB]
    vBC = [Point2f0(x) for x in BC]
    vCA = [Point2f0(x) for x in CA]
    poly!(scene, polyA, color = RGBf0(96/255,173/255,81/255))
    poly!(scene, polyB, color = RGBf0(170/255,121/255,193/255))
    poly!(scene, polyC, color = RGBf0(213/255,99/255,92/255))
    linesegments!(scene, vAB, color = RGBf0(102/255,130/225,223/255), linewidth = 5.0)
    linesegments!(scene, vBC, color = RGBf0(102/255,130/225,223/255), linewidth = 5.0)
    linesegments!(scene, vCA, color = RGBf0(102/255,130/225,223/255), linewidth = 5.0)
    scatter!(scene, vAB, color = RGBf0(64/255,99/225,216/255))
    scatter!(scene, vBC, color = RGBf0(64/255,99/225,216/255))
    scatter!(scene, vCA, color = RGBf0(64/255,99/225,216/255))
end

limits = FRect(-1, -2, 1, 1)
scene = Scene(resolution = (1000, 1000),
              show_axis = false, scale_plot = false,
              limits = limits)
drawjulialogo(scene)
