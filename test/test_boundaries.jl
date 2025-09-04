# test/test_boundaries.jl
using Test
using TriCo

# Half-open bins in r and μ; edge/near-edge behavior.
@testset "Boundary behavior (half-open bins)" begin
    rmin, rmax = 1.0, 2.0
    Nr, μmax, Nμ = 4, 1.0, 2

    # Three collinear points along x-axis:
    # distances: r12 ≈ rmin+ε, r13 ≈ rmax-ε, r23 ≈ (rmax - rmin) - small
    X = [0.0, rmin + 1e-9, rmax - 1e-9]
    Y = [0.0, 0.0, 0.0]
    Z = [0.0, 0.0, 0.0]

    # With μmax=1.0, μ=1 passes since we treat μ-bin via floor on sqrt(μ^2) ∈ [0,1)
    H = TriCo.count_triangles!(X,Y,Z; rmin=rmin, rmax=rmax, Nr=Nr, μmax=μmax, Nμ=Nμ, cellsize=rmax)
    @test sum(H.h) == 1
end

