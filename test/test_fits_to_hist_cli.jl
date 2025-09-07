#!/usr/bin/env julia
# test/test_fits_to_hist_cli.jl
using Test

@testset "fits_to_hist.jl mixed ABC runs and writes output" begin
    # Resolve paths relative to the package root (.. from test/)
    root   = normpath(joinpath(@__DIR__, ".."))
    script = joinpath(root, "scripts", "fits_to_hist.jl")

    a = joinpath(root, "data", "galaxies_A.fits")
    b = joinpath(root, "data", "galaxies_B.fits")
    c = joinpath(root, "data", "galaxies_C.fits")
    out = joinpath(root, "data", "galaxies.npz")

    # Clean up any previous output
    isfile(out) && rm(out; force=true)

    # Build the command
    cmd = `julia --project=. $script --fitsA $a --fitsB $b --fitsC $c \
                 --pattern ABC --rmin 5 --rmax 60 --Nr 55 \
                 --mumax 0.9 --Nmu 2 --cellsize 60 --out $out`

    # Run: throws if exit code ≠ 0
    run(cmd)

    # Assertions: file created and non-empty
    @test isfile(out)
    @test filesize(out) > 0

    # Optional: clean up after success (uncomment if desired)
    # rm(out; force=true)
end

