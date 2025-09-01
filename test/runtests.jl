import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))  # activate TriCo's project

using Test
using TriCo

include("test_geometry.jl")

