using Test, Random
using TriCo

@testset "Periodic basics (z is LOS; min-image)" begin
    Random.seed!(2)
    N = 400
    Lx, Ly, Lz = 100.0, 120.0, 80.0
    X, Y, Z = Lx .* rand(N), Ly .* rand(N), Lz .* rand(N)

    rmin, rmax, Nr = 5.0, 30.0, 5
    μmax, Nμ = 0.9, 4
    cellsize = rmax

    H = count_triangles_periodic!(X, Y, Z;
        Lx=Lx, Ly=Ly, Lz=Lz,
        rmin=rmin, rmax=rmax, Nr=Nr,
        μmax=μmax, Nμ=Nμ, cellsize=cellsize)

    @test size(H.h) == (Nr, Nr, Nμ, Nμ)
    @test sum(H.h) > 0
    @test all(>=(0), H.h)

    # Periodic shift invariance (wrap after integer box shifts)
    X2 = (X .+ Lx) .% Lx
    Y2 = (Y .+ 0.0) .% Ly
    Z2 = (Z .+ 2Lz) .% Lz
    Hshift = count_triangles_periodic!(X2, Y2, Z2;
        Lx=Lx, Ly=Ly, Lz=Lz,
        rmin=rmin, rmax=rmax, Nr=Nr,
        μmax=μmax, Nμ=Nμ, cellsize=cellsize)

    @test sum(Hshift.h) == sum(H.h)
end

