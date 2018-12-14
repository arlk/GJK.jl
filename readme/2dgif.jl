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

centerA = [ 0. √2.]
centerB = [ 1.  0.]
centerC = [-1.  0.]
polyA = randomcvx(15, centerA)
polyB = randomcvx(15, centerB)
polyC = randomcvx(15, centerC)
AB = closest_points(polyA', polyB', SVector{2}(centerA - centerB))
#  poly(polyA, color = RGBf0(51/255,151/255,31/255))
#  poly!(polyB, color = RGBf0(147/255,84/255,180/255))
#  poly!(polyC, color = RGBf0(204/255,51/255,26/255))
lines!([AB[1][1],AB[2][1]], [AB[1][2],AB[2][2]])
