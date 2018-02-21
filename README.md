# GJK

[![Build Status](https://travis-ci.com/arlk/GJK.jl.svg?token=r5fpy4W1YP6bTrztNxWQ&branch=master)](https://travis-ci.com/arlk/GJK.jl)

GJK.jl is a pure Julia package, designed to provide a full set of tools to determine several different collision detection queries between convex polytopes in 2/3 dimensions. Note, that this is a work in progress and not everything has been implemented yet.

The package allows you to work with polytopes defined as an array of vertices (Using StaticArrays will considerably speed up the query evaluations), for example:
```julia
polyA = @SMatrix rand(100, 2)
polyB = @SMatrix rand(100, 2) + 1.5
```

The collision detection query can be performed simply by providing the polytope information and an initial search direction as:
```julia
ret = gjk(polyA, polyB, SVector(0.0, 1.0))
```

[![codecov](https://codecov.io/gh/arlk/GJK.jl/branch/master/graph/badge.svg?token=f5dBcv0pe3)](https://codecov.io/gh/arlk/GJK.jl)

![](https://github.com/arlk/GJK.jl/raw/master/readme/collision2d.gif)

![](https://github.com/arlk/GJK.jl/raw/master/readme/collision3d.png)
