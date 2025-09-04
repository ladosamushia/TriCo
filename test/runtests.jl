using Test
using TriCo

@info "Running TriCo tests…"
include("test_nonperiodic.jl")   # keep your existing basics, if you want
include("test_periodic.jl")      # keep your existing basics, if you want
include("test_bruteforce.jl")    # the new reference checks

