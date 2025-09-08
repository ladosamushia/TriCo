# src/pairs_mixed.jl (updated to use pairs_utils)
# Cross-pair counting between two catalogs with (r, μ) binning and μ_max.
# Non-periodic uses true LOS; periodic uses z-LOS with min-image.

using Base.Threads
using .PairsUtils: HistRMu, mu_true_los_numden2, mu_z_los_numden2, mu_from_numden2,
                   bin_index_inv, bin_pair_pre!

# Expect build_grid_* and cell_points from triangles.jl.

if !isdefined(@__MODULE__, :cell_points)
    @inline function cell_points(starts::Vector{Int}, counts::Vector{Int}, pts::Vector{Int}, cid::Int)
        s = starts[cid]; c = counts[cid]
        return @view pts[s : s + c - 1]
    end
end

function count_pairs_cross_grid!(X1::AbstractVector{<:Real},
                                 Y1::AbstractVector{<:Real},
                                 Z1::AbstractVector{<:Real},
                                 X2::AbstractVector{<:Real},
                                 Y2::AbstractVector{<:Real},
                                 Z2::AbstractVector{<:Real};
                                 rmin::Real, rmax::Real, Nr::Integer,
                                 μmax::Real, Nμ::Integer=2,
                                 cellsize::Real=rmax)

    N1 = length(X1); @assert length(Y1) == N1 && length(Z1) == N1
    N2 = length(X2); @assert length(Y2) == N2 && length(Z2) == N2

    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf * μmaxf
    NrI = Int(Nr); NμI = Int(Nμ)

    invΔr = NrI / (rmaxf - rminf)
    invΔμ = NμI / μmaxf
    r2min = rminf * rminf
    r2max = rmaxf * rmaxf

    @assert isdefined(@__MODULE__, :build_grid_nonperiodic) "build_grid_nonperiodic not found."
    G = build_grid_nonperiodic(X2, Y2, Z2, cellsize)  # grid over catalog 2

    Hs = [zeros(Int, NrI, NμI) for _ in 1:Threads.nthreads()]

    @inbounds Threads.@threads for i in 1:N1
        tid = Threads.threadid()
        Hh  = Hs[tid]

        xi = Float64(X1[i]); yi = Float64(Y1[i]); zi = Float64(Z1[i])

        ix = Int(floor((xi - G.x0) * G.invcs))
        iy = Int(floor((yi - G.y0) * G.invcs))
        iz = Int(floor((zi - G.z0) * G.invcs))

        x0 = max(ix - 1, 0); x1 = min(ix + 1, G.nx - 1)
        y0 = max(iy - 1, 0); y1 = min(iy + 1, G.ny - 1)
        z0 = max(iz - 1, 0); z1 = min(iz + 1, G.nz - 1)

        for cx in x0:x1, cy in y0:y1, cz in z0:z1
            cid = ((cx * G.ny) + cy) * G.nz + cz + 1
            for j in cell_points(G.starts, G.counts, G.pts, cid)
                xj = Float64(X2[j]); yj = Float64(Y2[j]); zj = Float64(Z2[j])
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

    H = HistRMu(NrI, NμI; T=Int)
    @inbounds for t in 1:length(Hs)
        H.h .+= Hs[t]
    end
    return H
end

function count_pairs_cross_periodic_grid!(X1::AbstractVector{<:Real},
                                          Y1::AbstractVector{<:Real},
                                          Z1::AbstractVector{<:Real},
                                          X2::AbstractVector{<:Real},
                                          Y2::AbstractVector{<:Real},
                                          Z2::AbstractVector{<:Real};
                                          Lx::Real, Ly::Real, Lz::Real,
                                          rmin::Real, rmax::Real, Nr::Integer,
                                          μmax::Real, Nμ::Integer=2,
                                          cellsize::Real=rmax)

    N1 = length(X1); @assert length(Y1) == N1 && length(Z1) == N1
    N2 = length(X2); @assert length(Y2) == N2 && length(Z2) == N2

    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf * μmaxf
    NrI = Int(Nr); NμI = Int(Nμ)

    invΔr = NrI / (rmaxf - rminf)
    invΔμ = NμI / μmaxf
    r2min = rminf * rminf
    r2max = rmaxf * rmaxf

    @assert isdefined(@__MODULE__, :build_grid_periodic) "build_grid_periodic not found."
    G = build_grid_periodic(X2, Y2, Z2; Lx=Lx, Ly=Ly, Lz=Lz, cellsize=cellsize)

    Lxf, Lyf, Lzf = G.Lx, G.Ly, G.Lz
    hx, hy, hz = G.hx, G.hy, G.hz
    csx = Lxf / G.nx; csy = Lyf / G.ny; csz = Lzf / G.nz
    invcsx = 1.0 / csx; invcsy = 1.0 / csy; invcsz = 1.0 / csz

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

    @inbounds Threads.@threads for i in 1:N1
        tid = Threads.threadid()
        Hh  = Hs[tid]

        xi = Float64(X1[i]); yi = Float64(Y1[i]); zi = Float64(Z1[i])

        ix = Int(floor(xi * invcsx)) % G.nx; (ix < 0) && (ix += G.nx)
        iy = Int(floor(yi * invcsy)) % G.ny; (iy < 0) && (iy += G.ny)
        iz = Int(floor(zi * invcsz)) % G.nz; (iz < 0) && (iz += G.nz)

        for dx in -1:1, dy in -1:1, dz in -1:1
            jx = wrap(ix + dx, G.nx)
            jy = wrap(iy + dy, G.ny)
            jz = wrap(iz + dz, G.nz)
            cid = ((jx * G.ny) + jy) * G.nz + jz + 1
            for j in cell_points(G.starts, G.counts, G.pts, cid)
                xj = Float64(X2[j]); yj = Float64(Y2[j]); zj = Float64(Z2[j])
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

    H = HistRMu(NrI, NμI; T=Int)
    @inbounds for t in 1:length(Hs)
        H.h .+= Hs[t]
    end
    return H
end

