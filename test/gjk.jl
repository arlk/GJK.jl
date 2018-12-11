using GeometryTypes: Point, Vec

import GJK.support
using GeometryTypes: HyperSphere, HyperRectangle, HyperCube
function GJK.support(rect::HyperRectangle{N, T}, dir::AbstractVector) where {N, T}
    normvec = [if x > 0 1.0 else -1.0 end for x in normalize(dir./rect.widths, Inf)]
    SVector{N}(rect.widths.*normvec/2.0 + rect.origin)
end

function GJK.support(cube::HyperCube{N, T}, dir::AbstractVector) where {N, T}
    normvec = [if x > 0 1.0 else -1.0 end for x in normalize(dir, Inf)]
    SVector{N}(cube.width.*normvec/2.0 + cube.origin)
end

function GJK.support(sphere::HyperSphere{N, T}, dir::AbstractVector) where {N, T}
    SVector{N}(sphere.center + sphere.r*normalize(dir, 2))
end

@testset "gjk (2D)" begin
    vec2 = SVector{2}([1.0, 0.0])
    @testset "array of vertices" begin
        rect2 = SMatrix{3,4}([0.0 1.0 1.0 0.0; 0.0 0.0 2.0 2.0; 1.0 1.0 1.0 1.0])
        transform2(θ,x,y) = SMatrix{2,3}([cos(θ) -sin(θ) x; sin(θ) cos(θ) y])

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,3,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == false
        @test all(ret[3] .≈ [1.0, 0.0])
        @test all(ret[4] .≈ [3.0, 0.0])

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,1,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,0.5,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,3,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == false
        @test all(ret[3] .≈ [1.0, 1.4142135623730951])
        @test all(ret[4] .≈ [1.585786437626905, 1.4142135623730954])

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,1,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,0.5,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == true

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,3,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == false
        @test all(ret[3] .≈ [1.914213562373095, 0.9142135623730953])
        @test all(ret[4] .≈ [2.0, 1.0])

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,1,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == true

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,0.5,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret[1] == true
    end # array of vertices

    @testset "sphere" begin
        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(3.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec2)
        @test ret[1] == false
        @test all(ret[3] .≈ [1.0, 0.0])
        @test all(ret[4] .≈ [2.0, 0.0])

        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(2.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec2)
        @test ret[1] == true

        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(1.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec2)
        @test ret[1] == true
    end # sphere

    @testset "rectangle" begin
        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(3.0, 0.0), Vec(1.0, 2.0))
        ret = gjk(rectangleA, rectangleB, vec2)
        @test ret[1] == false
        @test all(ret[3] .≈ [0.5, -1.0])
        @test all(ret[4] .≈ [2.5, -1.0])

        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(1.0, 0.0), Vec(1.0, 2.0))
        ret = gjk(rectangleA, rectangleB, vec2)
        @test ret[1] == true

        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(0.5, 0.0), Vec(1.0, 2.0))
        ret = gjk(rectangleA, rectangleB, vec2)
        @test ret[1] == true
    end # rectangle

    @testset "cube" begin
        cubeA = HyperCube(Vec(0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(3.0, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec2)
        @test ret[1] == false
        @test all(ret[3] .≈ [0.5, -0.5])
        @test all(ret[4] .≈ [2.5, -0.5])

        cubeA = HyperCube(Vec(0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(1.0, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec2)
        @test ret[1] == true

        cubeA = HyperCube(Vec(0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(0.5, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec2)
        @test ret[1] == true
    end # cube
end # gjk

@testset "gjk (3D)" begin
    vec3 = SVector{3}([1.0, 0.0, 0.0])
    @testset "array of vertices" begin
        # TODO (add several tetrahedron configurations)
    end # array of vertices

    @testset "sphere" begin
        sphereA = HyperSphere(Point(0.0, 0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(3.0, 3.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec3)
        @test ret[1] == false
        @test all(ret[3] .≈ [0.7071067841796789, 0.7071067781934159, 0.0])
        @test all(ret[4] .≈ [2.2928932158203210, 2.2928932218065840, 0.0])

        sphereA = HyperSphere(Point(0.0, 0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(2.0, 0.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec3)
        @test ret[1] == true

        sphereA = HyperSphere(Point(0.0, 0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(1.0, 0.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec3)
        @test ret[1] == true
    end # sphere

    @testset "rectangle" begin
        rectangleA = HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(1.0, 2.0, 3.0))
        rectangleB = HyperRectangle(Vec(3.0, 3.0, 0.0), Vec(1.0, 2.0, 3.0))
        ret = gjk(rectangleA, rectangleB, vec3)
        @test ret[1] == false
        @test all(ret[3] .≈ [0.5, 1.0, -1.5])
        @test all(ret[4] .≈ [2.5, 2.0, -1.5])

        rectangleA = HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(1.0, 2.0, 3.0))
        rectangleB = HyperRectangle(Vec(1.0, 2.0, 0.0), Vec(1.0, 2.0, 3.0))
        ret = gjk(rectangleA, rectangleB, vec3)
        @test ret[1] == true

        rectangleA = HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(1.0, 2.0, 3.0))
        rectangleB = HyperRectangle(Vec(0.5, 0.0, 0.0), Vec(1.0, 2.0, 3.0))
        ret = gjk(rectangleA, rectangleB, vec3)
        @test ret[1] == true
    end # rectangle

    @testset "cube" begin
        cubeA = HyperCube(Vec(0.0, 0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(3.0, 3.0, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec3)
        @test ret[1] == false
        @test all(ret[3] .≈ [0.5, 0.5, -0.5])
        @test all(ret[4] .≈ [2.5, 2.5, -0.5])

        cubeA = HyperCube(Vec(0.0, 0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(1.0, 1.0, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec3)
        @test ret[1] == true

        cubeA = HyperCube(Vec(0.0, 0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(0.5, 0.0, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec3)
        @test ret[1] == true
    end # cube
end # gjk
