using Test, Random
using TriCo

@testset "Non-periodic basics" begin
    Random.seed!(0)
    N = 80
    X, Y, Z = randn(N), randn(N), randn(N)

    H = count_triangles!(X, Y, Z; rmin=0.0, rmax=1e9, Nr=1, μmax=1.0, Nμ=1, cellsize=1e9)
    @test size(H.h) == (1, 1, 1, 1)
    @test sum(H.h) == binomial(N, 3)

    # Totals invariant under permuting inputs
    Hπ = count_triangles!(Y, Z, X; rmin=0.0, rmax=1e9, Nr=1, μmax=1.0, Nμ=1, cellsize=1e9)
    @test sum(Hπ.h) == sum(H.h)
end

@testset "Non-periodic finite bins" begin
    Random.seed!(1)
    N = 120
    X, Y, Z = randn(N), randn(N), randn(N)
    H = count_triangles!(X, Y, Z; rmin=0.1, rmax=2.5, Nr=6, μmax=0.9, Nμ=5, cellsize=2.5)

    @test ndims(H.h) == 4
    @test size(H.h) == (6, 6, 5, 5)

    # sanity: all counts nonnegative, nothing NaN
    @test all(>=(0), H.h)
    @test !any(isnan, Float64.(H.h))
end

