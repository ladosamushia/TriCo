# src/pairs.jl (updated to use pairs_utils)
# Pair counting with (r, μ) binning and a μ_max threshold.
# Non-periodic (true LOS) and periodic (z-LOS) grid-accelerated versions.

using Base.Threads
using .PairsUtils: HistRMu, mu_true_los_numden2, mu_z_los_numden2, mu_from_numden2,
                   bin_index_inv, bin_pair_pre!

# Expecting build_grid_nonperiodic / build_grid_periodic and cell_points from triangles.jl.

# Safe access helper (only if triangles.jl didn't define it)
if !isdefined(@__MODULE__, :cell_points)
    @inline function cell_points(starts::Vector{Int}, counts::Vector{Int}, pts::Vector{Int}, cid::Int)
        s = starts[cid]; c = counts[cid]
        return @view pts[s : s + c - 1]
    end
end

function count_pairs_grid!(X::AbstractVector{<:Real},
                           Y::AbstractVector{<:Real},
                           Z::AbstractVector{<:Real};
                           rmin::Real, rmax::Real, Nr::Integer,
                           μmax::Real, Nμ::Integer=2,
                           cellsize::Real=rmax)

    N = length(X)
    @assert length(Y) == N && length(Z) == N

    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf * μmaxf
    NrI = Int(Nr); NμI = Int(Nμ)

    invΔr = NrI / (rmaxf - rminf)
    invΔμ = NμI / μmaxf
    r2min = rminf * rminf
    r2max = rmaxf * rmaxf

    @assert isdefined(@__MODULE__, :build_grid_nonperiodic) "build_grid_nonperiodic not found."
    G = build_grid_nonperiodic(X, Y, Z, cellsize)

    Hs = [zeros(Int, NrI, NμI) for _ in 1:Threads.nthreads()]

    @inbounds Threads.@threads for ci in 1:(G.nx*G.ny*G.nz)
        tid = Threads.threadid()
        Hh  = Hs[tid]

        ix = (ci - 1) ÷ (G.ny * G.nz)
        rem1 = (ci - 1) % (G.ny * G.nz)
        iy = rem1 ÷ G.nz
        iz = rem1 % G.nz

        pts_i = cell_points(G.starts, G.counts, G.pts, ci)
        length(pts_i) == 0 && continue

        x0 = max(ix - 1, 0); x1 = min(ix + 1, G.nx - 1)
        y0 = max(iy - 1, 0); y1 = min(iy + 1, G.ny - 1)
        z0 = max(iz - 1, 0); z1 = min(iz + 1, G.nz - 1)

        for idx_i in pts_i
            xi = Float64(X[idx_i]); yi = Float64(Y[idx_i]); zi = Float64(Z[idx_i])

            @inbounds for cx in x0:x1, cy in y0:y1, cz in z0:z1
                cid = ((cx * G.ny) + cy) * G.nz + cz + 1
                for idx_j in cell_points(G.starts, G.counts, G.pts, cid)
                    j = idx_j
                    (j <= idx_i) && continue
                    xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
                    dx = xj - xi; dy = yj - yi; dz = zj - zi
                    r2 = dx*dx + dy*dy + dz*dz
                    (r2 < r2min || r2 >= r2max) && continue
                    num2, den2 = mu_true_los_numden2(xi, yi, zi, xj, yj, zj)
                    !(den2 == 0.0 || num2 < μmax2 * den2) && continue
                    μ = mu_from_numden2(num2, den2)
                    r = sqrt(r2)
                    bin_pair_pre!(Hh, r, μ, rminf, invΔr, NrI, invΔμ, NμI)
                end
            end
        end
    end

    H = HistRMu(NrI, NμI; T=Int)
    @inbounds for t in 1:length(Hs)
        H.h .+= Hs[t]
    end
    return H
end

function count_pairs_periodic_grid!(X::AbstractVector{<:Real},
                                    Y::AbstractVector{<:Real},
                                    Z::AbstractVector{<:Real};
                                    Lx::Real, Ly::Real, Lz::Real,
                                    rmin::Real, rmax::Real, Nr::Integer,
                                    μmax::Real, Nμ::Integer=2,
                                    cellsize::Real=rmax)

    N = length(X)
    @assert length(Y) == N && length(Z) == N

    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf * μmaxf
    NrI = Int(Nr); NμI = Int(Nμ)

    invΔr = NrI / (rmaxf - rminf)
    invΔμ = NμI / μmaxf
    r2min = rminf * rminf
    r2max = rmaxf * rmaxf

    @assert isdefined(@__MODULE__, :build_grid_periodic) "build_grid_periodic not found."
    G = build_grid_periodic(X, Y, Z; Lx=Lx, Ly=Ly, Lz=Lz, cellsize=cellsize)
    Lxf, Lyf, Lzf = G.Lx, G.Ly, G.Lz
    hx, hy, hz = G.hx, G.hy, G.hz

    Hs = [zeros(Int, NrI, NμI) for _ in 1:Threads.nthreads()]

    @inline function minimg_d(Δ::Float64, halfL::Float64, L::Float64)
        if Δ >  halfL
            Δ - L
        elseif Δ <= -halfL
            Δ + L
        else
            Δ
        end
    end

    @inline function wrap(a::Int, n::Int)
        v = a % n
        (v < 0) && (v += n)
        return v
    end

    @inbounds Threads.@threads for ci in 1:(G.nx*G.ny*G.nz)
        tid = Threads.threadid()
        Hh  = Hs[tid]

        ix = (ci - 1) ÷ (G.ny * G.nz)
        rem1 = (ci - 1) % (G.ny * G.nz)
        iy = rem1 ÷ G.nz
        iz = rem1 % G.nz

        pts_i = cell_points(G.starts, G.counts, G.pts, ci)
        length(pts_i) == 0 && continue

        for idx_i in pts_i
            xi = Float64(X[idx_i]); yi = Float64(Y[idx_i]); zi = Float64(Z[idx_i])

            for dx in -1:1, dy in -1:1, dz in -1:1
                jx = wrap(ix + dx, G.nx)
                jy = wrap(iy + dy, G.ny)
                jz = wrap(iz + dz, G.nz)
                cid = ((jx * G.ny) + jy) * G.nz + jz + 1
                for idx_j in cell_points(G.starts, G.counts, G.pts, cid)
                    j = idx_j
                    (j <= idx_i) && continue
                    xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
                    dxij = minimg_d(xj - xi, hx, Lxf)
                    dyij = minimg_d(yj - yi, hy, Lyf)
                    dzij = minimg_d(zj - zi, hz, Lzf)
                    r2 = dxij*dxij + dyij*dyij + dzij*dzij
                    (r2 < r2min || r2 >= r2max) && continue
                    num2, den2 = mu_z_los_numden2(dxij, dyij, dzij)
                    !(den2 == 0.0 || num2 < μmax2 * den2) && continue
                    μ = mu_from_numden2(num2, den2)
                    r = sqrt(r2)
                    bin_pair_pre!(Hh, r, μ, rminf, invΔr, NrI, invΔμ, NμI)
                end
            end
        end
    end

    H = HistRMu(NrI, NμI; T=Int)
    @inbounds for t in 1:length(Hs)
        H.h .+= Hs[t]
    end
    return H
end

