using GeometryTypes: Point, Vec
using GeometryTypes: HyperSphere, HyperRectangle, HyperCube
using GJK: support

@testset "support functions" begin
    @testset "array of vertices" begin
        vertices = [3.0 5.0 7.0; 13.0 15.0 11.0]
        dir = [1.0, 2.0]
        @test all(support(vertices, dir) .≈ [5.0, 15.0])
        dir = [-1.0, 2.0]
        @test all(support(vertices, dir) .≈ [5.0, 15.0])
        dir = [-1.0, -2.0]
        @test all(support(vertices, dir) .≈ [3.0, 13.0])
        dir = [1.0, -2.0]
        @test all(support(vertices, dir) .≈ [7.0, 11.0])

        vertices = [3.0 5.0 7.0 5.0; 13.0 15.0 11.0 13.0; 1.0 2.0 3.0 4.0]
        dir = [1.0, 1.0, 2.0]
        @test all(support(vertices, dir) .≈ [5.0, 13.0, 4.0])
        dir = [1.0, 1.0, -2.0]
        @test all(support(vertices, dir) .≈ [5.0, 15.0, 2.0])
        dir = [1.0, -1.0, -2.0]
        @test all(support(vertices, dir) .≈ [7.0, 11.0, 3.0])
        dir = [-1.0, -1.0, -2.0]
        @test all(support(vertices, dir) .≈ [3.0, 13.0, 1.0])
    end # array of vertices

    @testset "sphere" begin
        sphere = HyperSphere(Point(1.0, 2.0), 10.0)
        dir = [1.0, 2.0]
        @test all(support(sphere, dir) .≈ [5.47213595499958, 10.94427190999916])
        dir = [-1.0, 2.0]
        @test all(support(sphere, dir) .≈ [-3.4721359549995796, 10.94427190999916])
        dir = [-1.0, -2.0]
        @test all(support(sphere, dir) .≈ [-3.4721359549995796, -6.944271909999159])
        dir = [1.0, -2.0]
        @test all(support(sphere, dir) .≈ [5.47213595499958, -6.944271909999159])

        sphere = HyperSphere(Point(1.0, 2.0, 3.0), 10.0)
        dir = [1.0, 1.0, 2.0]
        @test all(support(sphere, dir) .≈ [5.08248290463863, 6.08248290463863, 11.16496580927726])
        dir = [1.0, 1.0, -2.0]
        @test all(support(sphere, dir) .≈ [5.08248290463863, 6.08248290463863, -5.16496580927726])
        dir = [1.0, -1.0, -2.0]
        @test all(support(sphere, dir) .≈ [5.08248290463863, -2.08248290463863, -5.16496580927726])
        dir = [-1.0, -1.0, -2.0]
        @test all(support(sphere, dir) .≈ [-3.08248290463863, -2.08248290463863, -5.16496580927726])
    end # sphere

    @testset "rectangle" begin
        rectangle = HyperRectangle(Vec(1.0, 2.0), Vec(10.0, 20.0))
        dir = [1.0, 2.0]
        @test all(support(rectangle, dir) .≈ [6.0, 12.0])
        dir = [-1.0, 2.0]
        @test all(support(rectangle, dir) .≈ [-4.0, 12.0])
        dir = [-1.0, -2.0]
        @test all(support(rectangle, dir) .≈ [-4.0, -8.0])
        dir = [1.0, -2.0]
        @test all(support(rectangle, dir) .≈ [6.0, -8.0])

        rectangle = HyperRectangle(Vec(1.0, 2.0, 3.0), Vec(10.0, 20.0, 30.0))
        dir = [1.0, 1.0, 2.0]
        @test all(support(rectangle, dir) .≈ [6.0, 12.0, 18.0])
        dir = [1.0, 1.0, -2.0]
        @test all(support(rectangle, dir) .≈ [6.0, 12.0, -12.0])
        dir = [1.0, -1.0, -2.0]
        @test all(support(rectangle, dir) .≈ [6.0, -8.0, -12.0])
        dir = [-1.0, -1.0, -2.0]
        @test all(support(rectangle, dir) .≈ [-4.0, -8.0, -12.0])
    end # rectangle

    @testset "cube" begin
        cube = HyperCube(Vec(1.0, 2.0), 10.0)
        dir = [1.0, 2.0]
        @test all(support(cube, dir) .≈ [6.0, 7.0])
        dir = [-1.0, 2.0]
        @test all(support(cube, dir) .≈ [-4.0, 7.0])
        dir = [-1.0, -2.0]
        @test all(support(cube, dir) .≈ [-4.0, -3.0])
        dir = [1.0, -2.0]
        @test all(support(cube, dir) .≈ [6.0, -3.0])

        cube = HyperCube(Vec(1.0, 2.0, 3.0), 10.0)
        dir = [1.0, 1.0, 2.0]
        @test all(support(cube, dir) .≈ [6.0, 7.0, 8.0])
        dir = [1.0, 1.0, -2.0]
        @test all(support(cube, dir) .≈ [6.0, 7.0, -2.0])
        dir = [1.0, -1.0, -2.0]
        @test all(support(cube, dir) .≈ [6.0, -3.0, -2.0])
        dir = [-1.0, -1.0, -2.0]
        @test all(support(cube, dir) .≈ [-4.0, -3.0, -2.0])
    end # cube
end # support functions
