using Test, Random
using TriCo

@testset "Coarser vs finer binning total counts" begin
    Random.seed!(4)
    N = 200
    X, Y, Z = randn(N), randn(N), randn(N)

    H1 = count_triangles!(X, Y, Z; rmin=1.0, rmax=10.0, Nr=3, μmax=0.8, Nμ=2, cellsize=10.0)
    H2 = count_triangles!(X, Y, Z; rmin=1.0, rmax=10.0, Nr=6, μmax=0.8, Nμ=4, cellsize=10.0)

    @test sum(H1.h) == sum(H2.h)  # totals must match regardless of binning resolution
end

