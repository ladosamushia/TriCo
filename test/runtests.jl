using Test
using TriCo

@info "Running TriCo tests…"
include("bruteforce_compare.jl")
include("test_mixed.jl")
include("test_random_cube.jl")
include("test_fits_to_hist_cli.jl")

