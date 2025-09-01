# test/test_geometry.jl
using Test
using DeltaTri   # your package

@testset "Geometry tests" begin

    @testset "Pairwise distances periodic" begin
        # simple 3D points
        r1 = (0.0, 0.0, 0.0)
        r2 = (1.0, 0.0, 0.0)
        r3 = (0.0, 1.0, 0.0)

        # call your geometry function
        r12_sq, r23_sq, r31_sq, μ12_sq, μ23_sq, μ31_sq =
            triplet_geometry_periodic(r1..., r2..., r3...)

        # check squared distances
        @test r12_sq ≈ 1.0
        @test r23_sq ≈ 2.0    # between (1,0,0) and (0,1,0)
        @test r31_sq ≈ 1.0
    end

    @testset "Cosines periodic" begin
        # trivial test: symmetric points
        r1 = (1.0, 0.0, 0.0)
        r2 = (0.0, 1.0, 0.0)
        r3 = (0.0, 0.0, 1.0)

        _, _, _, μ12_sq, μ23_sq, μ31_sq =
            triplet_geometry_periodic(r1..., r2..., r3...)

        # Computed by hand
        @test μ12_sq ≈ 1.0
        @test μ23_sq ≈ 0.5
        @test μ31_sq ≈ 0.5
    end

    @testset "Pairwise distances periodic wrapping" begin
        # simple 3D points
        r1 = (0.0, 0.0, 0.0)
        r2 = (2001.0, 0.0, 0.0)
        r3 = (0.0, 2001.0, 0.0)

        # call your geometry function
        r12_sq, r23_sq, r31_sq, μ12_sq, μ23_sq, μ31_sq =
            triplet_geometry_periodic(r1..., r2..., r3...)

        # check squared distances
        @test r12_sq ≈ 1.0
        @test r23_sq ≈ 2.0    # between (1,0,0) and (0,1,0)
        @test r31_sq ≈ 1.0
    end

    @testset "Cosines periodic wrapping" begin
        # trivial test: symmetric points
        r1 = (1.0, 0.0, 0.0)
        r2 = (0.0, 2001.0, 0.0)
        r3 = (0.0, 0.0, 2001.0)

        _, _, _, μ12_sq, μ23_sq, μ31_sq =
            triplet_geometry_periodic(r1..., r2..., r3...)

        # Computed by hand
        @test μ12_sq ≈ 1.0
        @test μ23_sq ≈ 0.5
        @test μ31_sq ≈ 0.5
    end

end

