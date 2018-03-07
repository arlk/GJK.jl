# GJK

[![Build Status](https://travis-ci.org/arlk/GJK.jl.svg?branch=master)](https://travis-ci.org/arlk/GJK.jl) [![codecov](https://codecov.io/gh/arlk/GJK.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/arlk/GJK.jl)

![](https://github.com/arlk/GJK.jl/raw/master/readme/collision2d.gif)

GJK.jl implements the Gilber-Johnson-Keerthi Algorithm from their seminal paper on fast collision detection. The following query types are planned for two convex polytopes:
 - Minimum distance computation (released)
 - Tolerance verification (coming up)

The following queries will be explored at a later date (hopefully by June):
 - Boolean collision detection (no update so far)
 - Continous collision detection (no update so far)

## Usage

The package allows you to work with polytopes defined as an array of vertices (Using StaticArrays will considerably speed up the query evaluations), for example:
```julia
polyA = @SMatrix rand(100, 2)
polyB = @SMatrix rand(100, 2) + 1.5
```

The collision detection query can be performed simply by providing the polytope information and an initial search direction as:
```julia
ret = gjk(polyA, polyB, SVector(0.0, 1.0))
```

If you want to use your custom convex objects, you can do so by extending the support function as:
```julia
import GJK.support
function GJK.support(obj::MyFancyShape{N, T}, dir::AbstractArray{T, 1}) where {T<:AbstractFloat}
  # do something
end
```

## Examples

Minimum distance computation in 2D:

![](https://github.com/arlk/GJK.jl/raw/master/readme/collision2d.png)

Minimum distance computation in 3D:

![](https://github.com/arlk/GJK.jl/raw/master/readme/collision3d.png)

## Related Packages

[EnhancedGJK.jl](https://github.com/rdeits/EnhancedGJK.jl)

## References

Gilbert, E. G., D. W. Johnson, and S. S. Keerthi. “A Fast Procedure for Computing the Distance between Complex Objects in Three-Dimensional Space.” IEEE Journal on Robotics and Automation 4, no. 2 (April 1988): 193–203. https://doi.org/10.1109/56.2083.
