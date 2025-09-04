using Test, Random
using TriCo

@testset "Periodic basics (z is LOS; min-image)" begin
    Random.seed!(2)
    N = 400
    Lx, Ly, Lz = 100.0, 120.0, 80.0
    X, Y, Z = Lx .* rand(N), Ly .* rand(N), Lz .* rand(N)

    H = count_triangles_periodic!(X, Y, Z;
        Lx=Lx, Ly=Ly, Lz=Lz,
        rmin=5.0, rmax=30.0, Nr=5,
        μmax=0.9, Nμ=4, cellsize=30.0)

    @test size(H.h) == (5, 5, 4, 4)
    @test sum(H.h) > 0

    # permutations shouldn't change totals
    Hπ = count_triangles_periodic!(Y, Z, X;
        Lx=Ly, Ly=Lz, Lz=Lx,   # keep coordinates consistent with axes
        rmin=5.0, rmax=30.0, Nr=5, μmax=0.9, Nμ=4, cellsize=30.0)
    @test sum(Hπ.h) == sum(H.h)
end

