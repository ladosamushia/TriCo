# test/test_bruteforce.jl

using Test
using Random
using TriCo

# ---------------- Helpers (same conventions as TriCo now) ----------------

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

@inline function mu_z_los(dx::Float64, dy::Float64, dz::Float64)::Float64
    r2 = dx*dx + dy*dy + dz*dz
    (r2 == 0.0) ? 0.0 : (abs(dz) / sqrt(r2))
end

@inline function two_shortest_with_mus(r12::Float64, r23::Float64, r13::Float64,
                                       μ12::Float64, μ23::Float64, μ13::Float64)
    a_r, a_μ = r12, μ12
    b_r, b_μ = r23, μ23
    c_r, c_μ = r13, μ13

    if b_r < a_r
        a_r, b_r = b_r, a_r
        a_μ, b_μ = b_μ, a_μ
    end
    if c_r < a_r
        a_r, c_r = c_r, a_r
        a_μ, c_μ = c_μ, a_μ
    end
    if c_r < b_r
        return a_r, c_r, a_μ, c_μ
    else
        return a_r, b_r, a_μ, b_μ
    end
end

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

@inline function bin_index(v::Float64, vmin::Float64, vmax::Float64, N::Int)::Int
    (v < vmin || v >= vmax) && return 0
    Δ = (vmax - vmin) / N
    i = Int(floor((v - vmin) / Δ)) + 1
    (i < 1 || i > N) && return 0
    return i
end

@inline function bin_triangle!(H::Array{Int,4},
                               r_a::Float64, r_b::Float64, μ_a::Float64, μ_b::Float64,
                               rmin::Float64, rmax::Float64, Nr::Int,
                               μmax::Float64, Nμ::Int)
    ir = bin_index(r_a, rmin, rmax, Nr); ir == 0 && return
    jr = bin_index(r_b, rmin, rmax, Nr); jr == 0 && return
    im = bin_index(μ_a, 0.0,  μmax, Nμ); im == 0 && return
    jm = bin_index(μ_b, 0.0,  μmax, Nμ); jm == 0 && return
    @inbounds H[ir, jr, im, jm] += 1
end

# ---------------- Brute-force reference builders ----------------

function ref_hist_nonperiodic(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64};
                              rmin::Float64, rmax::Float64, Nr::Int,
                              μmax::Float64, Nμ::Int)
    N = length(X)
    H = zeros(Int, Nr, Nr, Nμ, Nμ)
    @inbounds for i in 1:N-2
        xi = X[i]; yi = Y[i]; zi = Z[i]
        for j in i+1:N-1
            xj = X[j]; yj = Y[j]; zj = Z[j]

            dxij = xj - xi; dyij = yj - yi; dzij = zj - zi
            r_ij = sqrt(dxij*dxij + dyij*dyij + dzij*dzij)
            (r_ij < rmin || r_ij >= rmax) && continue
            μ_ij = mu_true_los(xi, yi, zi, xj, yj, zj)
            μ_ij >= μmax && continue

            for k in j+1:N
                xk = X[k]; yk = Y[k]; zk = Z[k]

                dxjk = xk - xj; dyjk = yk - yj; dzjk = zk - zj
                dxik = xi - xk; dyik = yi - yk; dzik = zi - zk

                r_jk = sqrt(dxjk*dxjk + dyjk*dyjk + dzjk*dzjk)
                (r_jk < rmin || r_jk >= rmax) && continue
                r_ik = sqrt(dxik*dxik + dyik*dyik + dzik*dzik)
                (r_ik < rmin || r_ik >= rmax) && continue

                μ_jk = mu_true_los(xj, yj, zj, xk, yk, zk)
                μ_jk >= μmax && continue
                μ_ik = mu_true_los(xi, yi, zi, xk, yk, zk)
                μ_ik >= μmax && continue

                r_a, r_b, μ_a, μ_b = two_shortest_with_mus(r_ij, r_jk, r_ik, μ_ij, μ_jk, μ_ik)
                bin_triangle!(H, r_a, r_b, μ_a, μ_b, rmin, rmax, Nr, μmax, Nμ)
            end
        end
    end
    return H
end

function ref_hist_periodic(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64};
                           Lx::Float64, Ly::Float64, Lz::Float64,
                           rmin::Float64, rmax::Float64, Nr::Int,
                           μmax::Float64, Nμ::Int)
    N = length(X)
    H = zeros(Int, Nr, Nr, Nμ, Nμ)
    @inbounds for i in 1:N-2
        xi = X[i]; yi = Y[i]; zi = Z[i]
        for j in i+1:N-1
            xj = X[j]; yj = Y[j]; zj = Z[j]

            dxij = minimg(xj - xi, Lx); dyij = minimg(yj - yi, Ly); dzij = minimg(zj - zi, Lz)
            r_ij = sqrt(dxij*dxij + dyij*dyij + dzij*dzij)
            (r_ij < rmin || r_ij >= rmax) && continue
            μ_ij = mu_z_los(dxij, dyij, dzij)
            μ_ij >= μmax && continue

            for k in j+1:N
                xk = X[k]; yk = Y[k]; zk = Z[k]

                dxjk = minimg(xk - xj, Lx); dyjk = minimg(yk - yj, Ly); dzjk = minimg(zk - zj, Lz)
                dxik = minimg(xi - xk, Lx); dyik = minimg(yi - yk, Ly); dzik = minimg(zi - zk, Lz)

                r_jk = sqrt(dxjk*dxjk + dyjk*dyjk + dzjk*dzjk)
                (r_jk < rmin || r_jk >= rmax) && continue
                r_ik = sqrt(dxik*dxik + dyik*dyik + dzik*dzik)
                (r_ik < rmin || r_ik >= rmax) && continue

                μ_jk = mu_z_los(dxjk, dyjk, dzjk)
                μ_jk >= μmax && continue
                μ_ik = mu_z_los(dxik, dyik, dzik)
                μ_ik >= μmax && continue

                r_a, r_b, μ_a, μ_b = two_shortest_with_mus(r_ij, r_jk, r_ik, μ_ij, μ_jk, μ_ik)
                bin_triangle!(H, r_a, r_b, μ_a, μ_b, rmin, rmax, Nr, μmax, Nμ)
            end
        end
    end
    return H
end

# ---------------- Tiny diff reporter (only if equality fails) ----------------

function report_first_diffs(H_fast::Array{Int,4}, H_ref::Array{Int,4}; limit::Int=24)
    @assert size(H_fast) == size(H_ref)
    Nr, _, Nμa, Nμb = size(H_fast)
    diffs = Int[]
    for ir in 1:Nr, jr in 1:Nr, im in 1:Nμa, jm in 1:Nμb
        if H_fast[ir,jr,im,jm] != H_ref[ir,jr,im,jm]
            push!(diffs, 1)
            if length(diffs) <= limit
                @info "diff" ir=ir jr=jr im=im jm=jm fast=H_fast[ir,jr,im,jm] ref=H_ref[ir,jr,im,jm] Δ=(H_fast[ir,jr,im,jm]-H_ref[ir,jr,im,jm])
            end
        end
    end
    println("Nonzero differing bins: ", length(diffs))
end

# ---------------- Tests ----------------

@testset "Brute-force agreement: non-periodic (true LOS, linear μ, two-shortest)" begin
    Random.seed!(1234)
    N = 30
    X = randn(N); Y = randn(N); Z = randn(N)

    rmin, rmax, Nr = 0.25, 1.8, 5
    μmax, Nμ       = 0.9, 4
    cellsize       = rmax  # unused here, for API parity

    H_fast = count_triangles!(X, Y, Z; rmin=rmin, rmax=rmax, Nr=Nr,
                              μmax=μmax, Nμ=Nμ, cellsize=cellsize).h
    H_ref  = ref_hist_nonperiodic(Float64.(X), Float64.(Y), Float64.(Z);
                                  rmin=rmin, rmax=rmax, Nr=Nr, μmax=μmax, Nμ=Nμ)

    @test sum(H_fast) == sum(H_ref) ||
          (report_first_diffs(H_fast, H_ref); false)
    @test H_fast == H_ref ||
          (report_first_diffs(H_fast, H_ref); false)
end

@testset "Brute-force agreement: periodic (min-image, z-LOS, linear μ, two-shortest)" begin
    Random.seed!(2025)
    N = 32
    Lx, Ly, Lz = 20.0, 22.0, 18.0
    X = Lx .* rand(N); Y = Ly .* rand(N); Z = Lz .* rand(N)

    rmin, rmax, Nr = 0.6, 4.4, 6
    μmax, Nμ       = 0.85, 5
    cellsize       = rmax

    H_fast = count_triangles_periodic!(X, Y, Z; Lx=Lx, Ly=Ly, Lz=Lz,
                                       rmin=rmin, rmax=rmax, Nr=Nr,
                                       μmax=μmax, Nμ=Nμ, cellsize=cellsize).h
    H_ref  = ref_hist_periodic(Float64.(X), Float64.(Y), Float64.(Z);
                               Lx=Lx, Ly=Ly, Lz=Lz,
                               rmin=rmin, rmax=rmax, Nr=Nr, μmax=μmax, Nμ=Nμ)

    @test sum(H_fast) == sum(H_ref) ||
          (report_first_diffs(H_fast, H_ref); false)
    @test H_fast == H_ref ||
          (report_first_diffs(H_fast, H_ref); false)
end

