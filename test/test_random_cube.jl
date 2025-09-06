#!/usr/bin/env julia
using Test
using Random
using TriCo
using Base.Threads

@testset "Random cube triangle counting" begin
    N        = 2000         # small so the test runs quickly
    L        = 2000.0
    rmin     = 5.0
    rmax     = 60.0
    Nr       = 20
    μmax     = 0.9
    Nμ       = 2
    CELL     = rmax
    SEED     = 12345
    PERIODIC = false        # toggle to true if you want to cover periodic branch

    Random.seed!(SEED)

    X = L .* rand(N)
    Y = L .* rand(N)
    Z = L .* rand(N)

    H = if PERIODIC
        count_triangles_periodic_grid!(X, Y, Z;
            Lx=L, Ly=L, Lz=L,
            rmin=rmin, rmax=rmax, Nr=Nr,
            μmax=μmax, Nμ=Nμ, cellsize=CELL)
    else
        count_triangles_grid!(X, Y, Z;
            rmin=rmin, rmax=rmax, Nr=Nr,
            μmax=μmax, Nμ=Nμ, cellsize=CELL)
    end

    # Smoke-test conditions
    @test size(H.h) == (Nr, Nr, Nμ, Nμ)
    @test sum(H.h) > 0
end

