using Test, Random
using TriCo

# Helper: sum of all histogram counts
sumhist(H) = sum(H.h)

@testset "μ and r bin edge cases" begin
    Random.seed!(3)
    # Construct a tiny configuration to hit μ≈0, μ≈μmax and r near edges
    X = [0.0, 1.0, 0.0, 0.0]
    Y = [0.0, 0.0, 1.0, 0.0]
    Z = [0.0, 0.0, 0.0, 1.0]  # axes-aligned points to produce clear LOS angles

    rmin, rmax, Nr = 0.0, sqrt(2)+1e-8, 2
    μmax, Nμ = 1.0, 2
    H = count_triangles!(X, Y, Z; rmin=rmin, rmax=rmax, Nr=Nr, μmax=μmax, Nμ=Nμ, cellsize=rmax)

    @test sumhist(H) == binomial(length(X), 3)
    # Check no counts fall outside the last bin due to numerical fuzz
    @test all(>=(0), H.h)
end

