# src/pairs_utils.jl
# Shared utilities for pair counting: histogram container, μ helpers, and binning helpers.
# Intended to be included by src/pairs.jl and src/pairs_mixed.jl.

# =========================== Histogram ===========================
module PairsUtils

export HistRMu, HistRMu!, mu_true_los_numden2, mu_z_los_numden2, mu_from_numden2,
       bin_index_inv, bin_pair_pre!

mutable struct HistRMu{T<:Integer}
    h::Array{T,2}   # (Nr, Nμ)
end

@inline HistRMu(Nr::Integer, Nμ::Integer; T::Type{<:Integer}=Int) =
    HistRMu(zeros(T, Nr, Nμ))

# For symmetry with triangles style (not strictly needed)
@inline function HistRMu!(H::HistRMu{T}, Nr::Integer, Nμ::Integer) where {T<:Integer}
    resize!(H.h, Nr * Nμ)  # no-op pattern, kept for compatibility
    return H
end

# =========================== μ helpers ===========================

# μ = |cosθ| with true LOS (pair midpoint), return (num^2, den^2)
@inline function mu_true_los_numden2(xi::Float64, yi::Float64, zi::Float64,
                                     xj::Float64, yj::Float64, zj::Float64)
    sx = xj - xi;  sy = yj - yi;  sz = zj - zi
    lx = xi + xj;  ly = yi + yj;  lz = zi + zj
    s2 = sx*sx + sy*sy + sz*sz
    l2 = lx*lx + ly*ly + lz*lz
    num2 = (sx*lx + sy*ly + sz*lz)^2
    den2 = s2 * l2
    return num2, den2
end

# μ = |Δz| / r, return (num^2, den^2)
@inline function mu_z_los_numden2(dx::Float64, dy::Float64, dz::Float64)
    r2 = dx*dx + dy*dy + dz*dz
    return dz*dz, r2
end

@inline function mu_from_numden2(num2::Float64, den2::Float64)
    (den2 == 0.0) && return 0.0
    return sqrt(num2 / den2)
end

# =========================== Binning helpers =====================

@inline function bin_index_inv(v::Float64, vmin::Float64, invΔ::Float64, N::Int)::Int
    t = (v - vmin) * invΔ
    (t < 0.0 || t >= N) && return 0
    i = Int(floor(t)) + 1
    return (1 <= i <= N) ? i : 0
end

@inline function bin_pair_pre!(Hh::AbstractArray{<:Integer,2},
                               r::Float64, μ::Float64,
                               rmin::Float64, invΔr::Float64, Nr::Int,
                               invΔμ::Float64, Nμ::Int)
    ir = bin_index_inv(r, rmin, invΔr, Nr); ir == 0 && return
    im = bin_index_inv(μ, 0.0,  invΔμ, Nμ); im == 0 && return
    @inbounds Hh[ir, im] += 1
end

end # module

