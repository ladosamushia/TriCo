# src/triangles.jl

using Base.Threads

# 4D joint histogram over (r_a, r_b, μ_a, μ_b) where r_a ≤ r_b are the two
# shortest sides of the triangle, and μ_a, μ_b are the corresponding |cosθ|
# values attached to those sides.
mutable struct HistR12R23Mu12Mu13{T<:Integer}
    h::Array{T,4}   # size (Nr, Nr, Nμ, Nμ)
end

@inline HistR12R23Mu12Mu13(Nr::Integer, Nμ::Integer; T::Type{<:Integer}=Int) =
    HistR12R23Mu12Mu13(zeros(T, Nr, Nr, Nμ, Nμ))

# ---------------------- Helpers ----------------------

# μ for a pair using TRUE line-of-sight (pair midpoint), linear in μ.
@inline function mu_true_los(xi::Float64, yi::Float64, zi::Float64,
                             xj::Float64, yj::Float64, zj::Float64)::Float64
    sx = xj - xi;  sy = yj - yi;  sz = zj - zi
    lx = xi + xj;  ly = yi + yj;  lz = zi + zj
    s2 = sx*sx + sy*sy + sz*sz
    l2 = lx*lx + ly*ly + lz*lz
    if s2 == 0.0 || l2 == 0.0
        return 0.0
    end
    abs(sx*lx + sy*ly + sz*lz) / sqrt(s2 * l2)
end

# μ for a pair using FIXED z-LOS, linear in μ.
@inline function mu_z_los(dx::Float64, dy::Float64, dz::Float64)::Float64
    r2 = dx*dx + dy*dy + dz*dz
    (r2 == 0.0) ? 0.0 : (abs(dz) / sqrt(r2))
end

# Choose the two shortest sides and carry their μ’s alongside.
# Returns (r_a, r_b, μ_a, μ_b) with r_a ≤ r_b chosen from {r12,r23,r13}.
@inline function two_shortest_with_mus(r12::Float64, r23::Float64, r13::Float64,
                                       μ12::Float64, μ23::Float64, μ13::Float64)
    a_r, a_μ = r12, μ12
    b_r, b_μ = r23, μ23
    c_r, c_μ = r13, μ13

    # make a_r the smallest
    if b_r < a_r
        a_r, b_r = b_r, a_r
        a_μ, b_μ = b_μ, a_μ
    end
    if c_r < a_r
        a_r, c_r = c_r, a_r
        a_μ, c_μ = c_μ, a_μ
    end
    # choose second smallest between b and c
    if c_r < b_r
        return a_r, c_r, a_μ, c_μ
    else
        return a_r, b_r, a_μ, b_μ
    end
end

# Min-image on a single coordinate.
@inline function minimg(Δ::Float64, L::Float64)::Float64
    t = Δ - round(Δ / L) * L
    if t <= -L/2
        t + L
    elseif t > L/2
        t - L
    else
        t
    end
end

# Map scalar v in [vmin, vmax) into bin 1..N (else 0). Linear edges.
@inline function bin_index(v::Float64, vmin::Float64, vmax::Float64, N::Int)::Int
    (v < vmin || v >= vmax) && return 0
    Δ = (vmax - vmin) / N
    # floor gives 0..N-1, then +1 → 1..N
    i = Int(floor((v - vmin) / Δ)) + 1
    # guard against tiny roundoff
    (i < 1 || i > N) && return 0
    return i
end

# Push a single triangle into the 4D histogram with left-closed bins.
@inline function bin_triangle!(Hh::AbstractArray{<:Integer,4},
                               r_a::Float64, r_b::Float64, μ_a::Float64, μ_b::Float64,
                               rmin::Float64, rmax::Float64, Nr::Int,
                               μmax::Float64, Nμ::Int)
    ir = bin_index(r_a, rmin, rmax, Nr); ir == 0 && return
    jr = bin_index(r_b, rmin, rmax, Nr); jr == 0 && return
    im = bin_index(μ_a, 0.0,  μmax, Nμ); im == 0 && return
    jm = bin_index(μ_b, 0.0,  μmax, Nμ); jm == 0 && return
    @inbounds Hh[ir, jr, im, jm] += 1
end

# ---------------------- Non-periodic ----------------------

"""
    count_triangles!(X, Y, Z; rmin, rmax, Nr, μmax, Nμ=2, cellsize=rmax)

Count triangles where **all three pairs** satisfy:
- `rmin ≤ r_ij < rmax`
- `μ_ij < μmax` with μ = |cos θ|, **linear in μ**

For each triangle, identify the **two shortest sides** and bin their lengths
`(r_a, r_b)` and the corresponding μ’s `(μ_a, μ_b)` into a 4D histogram
with shape `(Nr, Nr, Nμ, Nμ)`. Bins are half-open `[min, max)`.

Returns a `HistR12R23Mu12Mu13`.
"""
function count_triangles!(X::AbstractVector{<:Real},
                          Y::AbstractVector{<:Real},
                          Z::AbstractVector{<:Real};
                          rmin::Real, rmax::Real, Nr::Integer,
                          μmax::Real, Nμ::Integer=2,
                          cellsize::Real=rmax)

    N = length(X)
    @assert length(Y) == N && length(Z) == N "X, Y, Z must have same length"
    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax)
    NrI   = Int(Nr)
    NμI   = Int(Nμ)

    H = HistR12R23Mu12Mu13(NrI, NμI; T=Int)

    # Simple O(N^3) loop; can be replaced by grid-accelerated traversal later.
    @inbounds for i in 1:N-2
        xi = Float64(X[i]); yi = Float64(Y[i]); zi = Float64(Z[i])
        for j in i+1:N-1
            xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])

            dxij = xj - xi; dyij = yj - yi; dzij = zj - zi
            r_ij = sqrt(dxij*dxij + dyij*dyij + dzij*dzij)
            (r_ij < rminf || r_ij >= rmaxf) && continue
            μ_ij = mu_true_los(xi, yi, zi, xj, yj, zj)
            μ_ij >= μmaxf && continue

            for k in j+1:N
                xk = Float64(X[k]); yk = Float64(Y[k]); zk = Float64(Z[k])

                dxjk = xk - xj; dyjk = yk - yj; dzjk = zk - zj
                dxik = xi - xk; dyik = yi - yk; dzik = zi - zk

                r_jk = sqrt(dxjk*dxjk + dyjk*dyjk + dzjk*dzjk)
                (r_jk < rminf || r_jk >= rmaxf) && continue
                r_ik = sqrt(dxik*dxik + dyik*dyik + dzik*dzik)
                (r_ik < rminf || r_ik >= rmaxf) && continue

                μ_jk = mu_true_los(xj, yj, zj, xk, yk, zk)
                μ_jk >= μmaxf && continue
                μ_ik = mu_true_los(xi, yi, zi, xk, yk, zk)
                μ_ik >= μmaxf && continue

                r_a, r_b, μ_a, μ_b = two_shortest_with_mus(r_ij, r_jk, r_ik, μ_ij, μ_jk, μ_ik)
                bin_triangle!(H.h, r_a, r_b, μ_a, μ_b, rminf, rmaxf, NrI, μmaxf, NμI)
            end
        end
    end

    return H
end

# ---------------------- Periodic ----------------------

"""
    count_triangles_periodic!(X, Y, Z;
        Lx, Ly, Lz, rmin, rmax, Nr, μmax, Nμ=2, cellsize=rmax)

Periodic-box variant with **minimum-image** separations and a fixed **z-LOS**.
Conventions match `count_triangles!`:

- Bin **linear in μ = |Δz|/r** (not μ²)
- Require **all three pairs** to satisfy `rmin ≤ r_ij < rmax` and `μ_ij < μmax`
- Bin the **two shortest sides** (and their μ’s) into `(Nr, Nr, Nμ, Nμ)` with half-open edges.

Returns a `HistR12R23Mu12Mu13`.
"""
function count_triangles_periodic!(X::AbstractVector{<:Real},
                                   Y::AbstractVector{<:Real},
                                   Z::AbstractVector{<:Real};
                                   Lx::Real, Ly::Real, Lz::Real,
                                   rmin::Real, rmax::Real, Nr::Integer,
                                   μmax::Real, Nμ::Integer=2,
                                   cellsize::Real=rmax)

    N = length(X)
    @assert length(Y) == N && length(Z) == N "X, Y, Z must have same length"
    Lxf = Float64(Lx); Lyf = Float64(Ly); Lzf = Float64(Lz)
    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax)
    NrI   = Int(Nr)
    NμI   = Int(Nμ)

    H = HistR12R23Mu12Mu13(NrI, NμI; T=Int)

    @inbounds for i in 1:N-2
        xi = Float64(X[i]); yi = Float64(Y[i]); zi = Float64(Z[i])
        for j in i+1:N-1
            xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])

            dxij = minimg(xj - xi, Lxf)
            dyij = minimg(yj - yi, Lyf)
            dzij = minimg(zj - zi, Lzf)
            r_ij = sqrt(dxij*dxij + dyij*dyij + dzij*dzij)
            (r_ij < rminf || r_ij >= rmaxf) && continue
            μ_ij = mu_z_los(dxij, dyij, dzij)
            μ_ij >= μmaxf && continue

            for k in j+1:N
                xk = Float64(X[k]); yk = Float64(Y[k]); zk = Float64(Z[k])

                dxjk = minimg(xk - xj, Lxf)
                dyjk = minimg(yk - yj, Lyf)
                dzjk = minimg(zk - zj, Lzf)
                dxik = minimg(xi - xk, Lxf)
                dyik = minimg(yi - yk, Lyf)
                dzik = minimg(zi - zk, Lzf)

                r_jk = sqrt(dxjk*dxjk + dyjk*dyjk + dzjk*dzjk)
                (r_jk < rminf || r_jk >= rmaxf) && continue
                r_ik = sqrt(dxik*dxik + dyik*dyik + dzik*dzik)
                (r_ik < rminf || r_ik >= rmaxf) && continue

                μ_jk = mu_z_los(dxjk, dyjk, dzjk)
                μ_jk >= μmaxf && continue
                μ_ik = mu_z_los(dxik, dyik, dzik)
                μ_ik >= μmaxf && continue

                r_a, r_b, μ_a, μ_b = two_shortest_with_mus(r_ij, r_jk, r_ik, μ_ij, μ_jk, μ_ik)
                bin_triangle!(H.h, r_a, r_b, μ_a, μ_b, rminf, rmaxf, NrI, μmaxf, NμI)
            end
        end
    end

    return H
end

