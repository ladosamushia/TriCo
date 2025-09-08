#!/usr/bin/env julia
using Test
using Random
using TriCo
using Base.Threads

# -------------------- Smoke test (fast) --------------------
@testset "Random cube pair counting (smoke)" begin
    N        = 2000         # small so the test runs quickly
    L        = 2000.0
    rmin     = 5.0
    rmax     = 60.0
    Nr       = 20
    μmax     = 0.9
    Nμ       = 2
    CELL     = rmax
    SEED     = 12345
    PERIODIC = false        # toggle to true if you want to cover periodic branch

    Random.seed!(SEED)

    X = L .* rand(N)
    Y = L .* rand(N)
    Z = L .* rand(N)

    H = if PERIODIC
        TriCo.count_pairs_periodic_grid!(X, Y, Z;
            Lx=L, Ly=L, Lz=L,
            rmin=rmin, rmax=rmax, Nr=Nr,
            μmax=μmax, Nμ=Nμ, cellsize=CELL)
    else
        TriCo.count_pairs_grid!(X, Y, Z;
            rmin=rmin, rmax=rmax, Nr=Nr,
            μmax=μmax, Nμ=Nμ, cellsize=CELL)
    end

    # Smoke-test conditions
    @test size(H.h) == (Nr, Nμ)
    @test sum(H.h) > 0
end

# -------------------- Brute-force agreement (non-periodic) --------------------
@testset "Brute-force agreement: non-periodic (true LOS)" begin
    N        = 600          # keep modest so O(N^2) runs fast
    L        = 500.0
    rmin     = 5.0
    rmax     = 60.0
    Nr       = 12
    μmax     = 0.9
    Nμ       = 3
    CELL     = rmax
    SEED     = 42

    Random.seed!(SEED)
    X = L .* rand(N); Y = L .* rand(N); Z = L .* rand(N)

    H = TriCo.count_pairs_grid!(X, Y, Z; rmin=rmin, rmax=rmax, Nr=Nr, μmax=μmax, Nμ=Nμ, cellsize=CELL)

    # ---- Brute force reference ----
    NrI = Int(Nr); NμI = Int(Nμ)
    rminf = Float64(rmin); rmaxf = Float64(rmax); invΔr = NrI / (rmaxf - rminf)
    μmaxf = Float64(μmax); invΔμ = NμI / μmaxf
    r2min = rminf * rminf; r2max = rmaxf * rmaxf
    μmax2 = μmaxf * μmaxf

    hist = zeros(Int, NrI, NμI)

    @inline function bin_index_inv(v::Float64, vmin::Float64, invΔ::Float64, N::Int)::Int
        t = (v - vmin) * invΔ
        (t < 0.0 || t >= N) && return 0
        i = Int(floor(t)) + 1
        return (1 <= i <= N) ? i : 0
    end

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

    @inbounds for i in 1:N-1, j in i+1:N
        xi = Float64(X[i]); yi = Float64(Y[i]); zi = Float64(Z[i])
        xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
        dx = xj - xi; dy = yj - yi; dz = zj - zi
        r2 = dx*dx + dy*dy + dz*dz
        (r2 < r2min || r2 >= r2max) && continue
        num2, den2 = mu_true_los_numden2(xi, yi, zi, xj, yj, zj)
        !(den2 == 0.0 || num2 < μmax2 * den2) && continue
        μ = den2 == 0.0 ? 0.0 : sqrt(num2/den2)
        r = sqrt(r2)
        ir = bin_index_inv(r, rminf, invΔr, NrI); ir == 0 && continue
        im = bin_index_inv(μ, 0.0,  invΔμ, NμI); im == 0 && continue
        @inbounds hist[ir, im] += 1
    end

    @test hist == H.h
end

# -------------------- Brute-force agreement (periodic) --------------------
@testset "Brute-force agreement: periodic (min-image, z-LOS)" begin
    N        = 500
    Lx       = 400.0
    Ly       = 500.0
    Lz       = 450.0
    rmin     = 5.0
    rmax     = 60.0
    Nr       = 12
    μmax     = 0.9
    Nμ       = 3
    CELL     = rmax
    SEED     = 7

    Random.seed!(SEED)
    X = Lx .* rand(N); Y = Ly .* rand(N); Z = Lz .* rand(N)

    H = TriCo.count_pairs_periodic_grid!(X, Y, Z; Lx=Lx, Ly=Ly, Lz=Lz,
                                         rmin=rmin, rmax=rmax, Nr=Nr, μmax=μmax, Nμ=Nμ, cellsize=CELL)

    NrI = Int(Nr); NμI = Int(Nμ)
    rminf = Float64(rmin); rmaxf = Float64(rmax); invΔr = NrI / (rmaxf - rminf)
    μmaxf = Float64(μmax); invΔμ = NμI / μmaxf
    r2min = rminf * rminf; r2max = rmaxf * rmaxf
    μmax2 = μmaxf * μmaxf

    hist = zeros(Int, NrI, NμI)

    Lxf, Lyf, Lzf = Float64(Lx), Float64(Ly), Float64(Lz)
    hx, hy, hz = 0.5Lxf, 0.5Lyf, 0.5Lzf

    @inline function minimg_d(Δ::Float64, halfL::Float64, L::Float64)
        if Δ >  halfL
            Δ - L
        elseif Δ <= -halfL
            Δ + L
        else
            Δ
        end
    end

    @inline function bin_index_inv(v::Float64, vmin::Float64, invΔ::Float64, N::Int)::Int
        t = (v - vmin) * invΔ
        (t < 0.0 || t >= N) && return 0
        i = Int(floor(t)) + 1
        return (1 <= i <= N) ? i : 0
    end

    @inline function mu_z_los_numden2(dx::Float64, dy::Float64, dz::Float64)
        r2 = dx*dx + dy*dy + dz*dz
        return dz*dz, r2
    end

    @inbounds for i in 1:N-1, j in i+1:N
        xi = Float64(X[i]); yi = Float64(Y[i]); zi = Float64(Z[i])
        xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
        dx = minimg_d(xj - xi, hx, Lxf)
        dy = minimg_d(yj - yi, hy, Lyf)
        dz = minimg_d(zj - zi, hz, Lzf)
        r2 = dx*dx + dy*dy + dz*dz
        (r2 < r2min || r2 >= r2max) && continue
        num2, den2 = mu_z_los_numden2(dx, dy, dz)
        !(den2 == 0.0 || num2 < μmax2 * den2) && continue
        μ = den2 == 0.0 ? 0.0 : sqrt(num2/den2)
        r = sqrt(r2)
        ir = bin_index_inv(r, rminf, invΔr, NrI); ir == 0 && continue
        im = bin_index_inv(μ, 0.0,  invΔμ, NμI); im == 0 && continue
        @inbounds hist[ir, im] += 1
    end

    @test hist == H.h
end

