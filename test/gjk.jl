using GeometryTypes: Point, Vec
using GeometryTypes: HyperSphere, HyperRectangle, HyperCube

@testset "gjk" begin
    vec2 = [1.0, 0.0]
    vec3 = [1.0, 0.0, 0.0]
    @testset "array of vertices" begin
        rect2 = [0.0 1.0 1.0 0.0; 0.0 0.0 2.0 2.0; 1.0 1.0 1.0 1.0]
        transform2(θ,x,y) = [cos(θ) -sin(θ) x; sin(θ) cos(θ) y]

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,3,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == false
        @test all(ret.contact .≈ [1.0 3.0; 0.0 0.0])

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,1,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,0.5,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,3,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == false
        @test all(ret.contact .≈ [1.0 1.585786437626905; 1.4142135623730951 1.4142135623730954])

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,1,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,0.5,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == true

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,3,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == false
        @test all(ret.contact .≈ [1.914213562373095 2.0; 0.9142135623730953 1.0])

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,1,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == true

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,0.5,0) * rect2
        ret = gjk(polyA, polyB, vec2)
        @test ret.status == true
        @test ret.collision == true
    end # array of vertices

    @testset "sphere" begin
        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(3.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec2)
        @test ret.status == true
        @test ret.collision == false
        @test all(ret.contact .≈ [1.0 2.0; 0.0 0.0])

        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(2.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec2)
        @test ret.status == true
        @test ret.collision == true

        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(1.0, 0.0), 1.0)
        ret = gjk(sphereA, sphereB, vec2)
        @test ret.status == true
        @test ret.collision == true
    end # sphere

    @testset "rectangle" begin
        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(3.0, 0.0), Vec(1.0, 2.0))
        ret = gjk(rectangleA, rectangleB, vec2)
        @test ret.status == true
        @test ret.collision == false
        @test all(ret.contact .≈ [0.5 2.5; -1.0 -1.0])

        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(1.0, 0.0), Vec(1.0, 2.0))
        ret = gjk(rectangleA, rectangleB, vec2)
        @test ret.status == true
        @test ret.collision == true

        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(0.5, 0.0), Vec(1.0, 2.0))
        ret = gjk(rectangleA, rectangleB, vec2)
        @test ret.status == true
        @test ret.collision == true
    end # rectangle

    @testset "cube" begin
        cubeA = HyperCube(Vec(0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(3.0, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec2)
        @test ret.status == true
        @test ret.collision == false
        @test all(ret.contact .≈ [0.5 2.5; -0.5 -0.5])

        cubeA = HyperCube(Vec(0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(1.0, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec2)
        @test ret.status == true
        @test ret.collision == true

        cubeA = HyperCube(Vec(0.0, 0.0), 1.0)
        cubeB = HyperCube(Vec(0.5, 0.0), 1.0)
        ret = gjk(cubeA, cubeB, vec2)
        @test ret.status == true
        @test ret.collision == true
    end # cube
end # gjk

