# GJK

[![Build Status](https://travis-ci.com/arlk/Collide.jl.svg?branch=master)](https://travis-ci.com/arlk/GJK.jl) [![codecov](https://codecov.io/gh/arlk/GJK.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/arlk/GJK.jl)

![](https://github.com/arlk/GJK.jl/raw/master/readme/collision2d.gif)

GJK.jl is a pure Julia package, designed to provide a full set of tools to determine several different collision detection queries between convex polytopes in 2/3 dimensions. Note, that this is a work in progress and not everything has been implemented yet.

# Roadmap
(TODO)

# Usage

The package allows you to work with polytopes defined as an array of vertices (Using StaticArrays will considerably speed up the query evaluations), for example:
```julia
polyA = @SMatrix rand(100, 2)
polyB = @SMatrix rand(100, 2) + 1.5
```

The collision detection query can be performed simply by providing the polytope information and an initial search direction as:
```julia
ret = gjk(polyA, polyB, SVector(0.0, 1.0))
```

# Examples

Minimum distance computation in 2D:

![](https://github.com/arlk/GJK.jl/raw/master/readme/collision2d.png)

Minimum distance computation in 3D:

![](https://github.com/arlk/GJK.jl/raw/master/readme/collision3d.png)
