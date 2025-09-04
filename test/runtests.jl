using Test
using TriCo

@info "Running TriCo tests…"
include("test_nonperiodic.jl")
include("test_periodic.jl")
include("test_edges.jl")
include("test_regression.jl")

