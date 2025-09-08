# test/bruteforce_compare_FIXED.jl
using Test
using Random
using LinearAlgebra

using TriCo

# ----- Helpers matching src/triangles.jl --------------------------------------------------------
@inline function mu_true_los(xi,yi,zi, xj,yj,zj)
    sx = xj - xi;  sy = yj - yi;  sz = zj - zi
    s2 = sx*sx + sy*sy + sz*sz
    lx = xi + xj;  ly = yi + yj;  lz = zi + zj
    l2 = lx*lx + ly*ly + lz*lz
    if s2 == 0.0 || l2 == 0.0
        return 0.0
    end
    return abs(sx*lx + sy*ly + sz*lz) / sqrt(s2*l2)
end

@inline function bin_index_right_open(v::Float64, vmin::Float64, invΔ::Float64, N::Int)::Int
    t = (v - vmin) * invΔ
    (t < 0.0 || t >= N) && return 0
    return Int(floor(t)) + 1
end

function make_edges(rmin::Real, rmax::Real, Nr::Integer, μmax::Real, Nμ::Integer;
                    rspacing::Symbol=:linear)
    r_edges = rspacing === :linear ? range(rmin, rmax; length=Nr+1) :
               rspacing === :log    ? exp.(range(log(max(rmin, eps())), log(rmax); length=Nr+1)) :
               error("Unsupported rspacing=$(rspacing).")
    μ_edges = range(0.0, μmax; length=Nμ+1)
    return collect(r_edges), collect(μ_edges)
end

# Return indices (ia, ib) of the two smallest of (a,b,c), tie-breaking deterministically.
@inline function two_smallest_indices(a::Float64, b::Float64, c::Float64)
    ia, ib, ic = 1, 2, 3
    ra, rb, rc = a, b, c
    if rb < ra || (rb == ra && ib < ia)
        ra, rb = rb, ra
        ia, ib = ib, ia
    end
    if rc < ra || (rc == ra && ic < ia)
        ra, rc = rc, ra
        ia, ic = ic, ia
    end
    if rc < rb || (rc == rb && ic < ib)
        return ia, ic   # (ra, rc)
    else
        return ia, ib   # (ra, rb)
    end
end

# Brute force that EXACTLY mirrors pruning and binning in count_triangles_grid!
function brute_mirror_nonperiodic_true(X::AbstractVector, Y::AbstractVector, Z::AbstractVector;
                                       rmin::Real, rmax::Real, Nr::Integer,
                                       μmax::Real, Nμ::Integer)
    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf*μmaxf
    invΔr = Int(Nr) / (rmaxf - rminf)
    invΔμ = Int(Nμ) / μmaxf
    r2min = rminf*rminf; r2max = rmaxf*rmaxf

    H = zeros(Int, Nr, Nr, Nμ, Nμ)

    Np = length(X)
    @inbounds for i in 1:Np-2
        xi = Float64(X[i]); yi = Float64(Y[i]); zi = Float64(Z[i])
        for j in i+1:Np-1
            xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
            dx = xj - xi; dy = yj - yi; dz = zj - zi
            r2_ij = dx*dx + dy*dy + dz*dz
            (r2_ij < r2min || r2_ij >= r2max) && continue
            num2 = (dx*(xi+xj) + dy*(yi+yj) + dz*(zi+zj))^2
            den2 = r2_ij * ((xi+xj)^2 + (yi+yj)^2 + (zi+zj)^2)
            !(den2 == 0.0 || num2 < μmax2*den2) && continue

            for k in j+1:Np
                xk = Float64(X[k]); yk = Float64(Y[k]); zk = Float64(Z[k])
                dxik = xk - xi; dyik = yk - yi; dzik = zk - zi
                r2_ik = dxik*dxik + dyik*dyik + dzik*dzik
                (r2_ik < r2min || r2_ik >= r2max) && continue
                num2 = (dxik*(xi+xk) + dyik*(yi+yk) + dzik*(zi+zk))^2
                den2 = r2_ik * ((xi+xk)^2 + (yi+yk)^2 + (zi+zk)^2)
                !(den2 == 0.0 || num2 < μmax2*den2) && continue

                dxjk = xk - xj; dyjk = yk - yj; dzjk = zk - zj
                r2_jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk
                (r2_jk < r2min || r2_jk >= r2max) && continue
                num2 = (dxjk*(xj+xk) + dyjk*(yj+yk) + dzjk*(zj+zk))^2
                den2 = r2_jk * ((xj+xk)^2 + (yj+yk)^2 + (zj+zk)^2)
                !(den2 == 0.0 || num2 < μmax2*den2) && continue

                # two shortest with μ paired correctly
                ia, ib = two_smallest_indices(r2_ij, r2_jk, r2_ik)
                μij = mu_true_los(xi,yi,zi, xj,yj,zj)
                μjk = mu_true_los(xj,yj,zj, xk,yk,zk)
                μik = mu_true_los(xi,yi,zi, xk,yk,zk)
                r2s = (r2_ij, r2_jk, r2_ik)
                μs  = (μij, μjk, μik)

                r_a = sqrt(r2s[ia]); r_b = sqrt(r2s[ib])
                μ_a = μs[ia];        μ_b = μs[ib]

                ir = bin_index_right_open(r_a, rminf, invΔr, Int(Nr)); ir==0 && continue
                jr = bin_index_right_open(r_b, rminf, invΔr, Int(Nr)); jr==0 && continue
                im = bin_index_right_open(μ_a, 0.0,   invΔμ, Int(Nμ)); im==0 && continue
                jm = bin_index_right_open(μ_b, 0.0,   invΔμ, Int(Nμ)); jm==0 && continue
                H[ir,jr,im,jm] += 1
            end
        end
    end
    return H
end

# ---------- Tests --------------------------------------------------------------------------------
@testset "Brute-force agreement: NON-PERIODIC, TRUE LOS, linear μ, two-shortest" begin
    Random.seed!(12345)
    N   = 42
    X   = randn(Float64, N); Y = randn(Float64, N); Z = randn(Float64, N)

    r_edges, μ_edges = make_edges(0.25, 2.7, 6, 0.9, 5; rspacing=:linear)
    rmin, rmax = first(r_edges), last(r_edges)
    Nr  = length(r_edges)-1
    μmax = last(μ_edges)
    Nμ  = length(μ_edges)-1

    H_brute = brute_mirror_nonperiodic_true(X,Y,Z; rmin=rmin, rmax=rmax, Nr=Nr, μmax=μmax, Nμ=Nμ)

    H_fast = count_triangles_grid!(X, Y, Z;
                                   rmin=rmin, rmax=rmax, Nr=Nr,
                                   μmax=μmax, Nμ=Nμ,
                                   cellsize=rmax).h

    @test size(H_fast) == size(H_brute)
    @test H_fast == H_brute
end

println("Non-periodic TRUE-LOS brute-force comparison passed.")

