using GeometryTypes: Point, Vec
using GeometryTypes: HyperSphere, HyperRectangle, HyperCube
using Collide: support

@testset "support mapping (2D)" begin
    @testset "array of vertices" begin
        vertices = [3.0 5.0 7.0; 13.0 15.0 11.0]
        dir = [1.0, 2.0]
        @test support(vertices, dir) .≈ [5.0, 15.0]
        dir = [-1.0, 2.0]
        @test support(vertices, dir) .≈ [5.0, 15.0]
        dir = [-1.0, -2.0]
        @test support(vertices, dir) .≈ [3.0, 13.0]
        dir = [1.0, -2.0]
        @test support(vertices, dir) .≈ [7.0, 11.0]
    end

    @testset "circle" begin
        circle = HyperSphere(Point(1.0, 2.0), 10.0)
        dir = [1.0, 2.0]
        @test support(circle, dir) .≈ [5.44721, 10.9443]
        dir = [-1.0, 2.0]
        @test support(circle, dir) .≈ [-3.47214, 10.9443]
        dir = [-1.0, -2.0]
        @test support(circle, dir) .≈ [-3.47214, -6.10557]
        dir = [1.0, -2.0]
        @test support(circle, dir) .≈ [5.44721, -6.10557]
    end

    @testset "rectangle" begin
        rectangle = HyperRectangle(Vec(1.0, 2.0), Vec(10.0, 20.0))
        dir = [1.0, 2.0]
        @test support(rectangle, dir) .≈ [6.0, 12.0]
        dir = [-1.0, 2.0]
        @test support(rectangle, dir) .≈ [-4.0, 12.0]
        dir = [-1.0, -2.0]
        @test support(rectangle, dir) .≈ [-4.0, -8.0]
        dir = [1.0, -2.0]
        @test support(rectangle, dir) .≈ [6.0, -8.0]
    end

    @testset "square" begin
        square = HyperCube(Vec(1.0, 2.0), 10.0)
        dir = [1.0, 2.0]
        @test support(square, dir) .≈ [3.5, 7.0]
        dir = [-1.0, 2.0]
        @test support(square, dir) .≈ [-1.5, 7.0]
        dir = [-1.0, -2.0]
        @test support(square, dir) .≈ [-1.5, -3.0]
        dir = [1.0, -2.0]
        @test support(square, dir) .≈ [3.5, -3.0]
    end
end
