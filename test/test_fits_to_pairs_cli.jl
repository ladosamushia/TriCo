#!/usr/bin/env julia
# test/test_fits_to_pairs_cli.jl
using Test

@testset "fits_to_pairs.jl CLI" begin
    # Resolve paths relative to the package root (.. from test/)
    root   = normpath(joinpath(@__DIR__, ".."))
    script = joinpath(root, "scripts", "fits_to_pairs.jl")

    a = joinpath(root, "data", "galaxies_A.fits")
    b = joinpath(root, "data", "galaxies_B.fits")

    @testset "cross AB runs and writes output" begin
        out_ab = joinpath(root, "data", "pairs_AB.npz")
        isfile(out_ab) && rm(out_ab; force=true)

        cmd = `julia --project=. $script --fitsA $a --fitsB $b \
                     --rmin 5 --rmax 60 --Nr 20 \
                     --mumax 0.9 --Nmu 2 --cellsize 60 --out $out_ab`

        run(cmd)

        @test isfile(out_ab)
        @test filesize(out_ab) > 0
        # rm(out_ab; force=true) # optional cleanup
    end

    @testset "single A runs and writes output" begin
        out_a = joinpath(root, "data", "pairs_A.npz")
        isfile(out_a) && rm(out_a; force=true)

        cmd = `julia --project=. $script --fits $a \
                     --rmin 5 --rmax 60 --Nr 20 \
                     --mumax 0.9 --Nmu 2 --cellsize 60 --out $out_a`

        run(cmd)

        @test isfile(out_a)
        @test filesize(out_a) > 0
        # rm(out_a; force=true) # optional cleanup
    end
end

