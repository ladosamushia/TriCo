# src/triangles.jl

using Base.Threads

using .PairsUtils: mu_true_los_numden2, mu_z_los_numden2, mu_from_numden2,
                   bin_index_inv

# =========================== Histogram ===========================

mutable struct HistR12R23Mu12Mu13{T<:Integer}
    h::Array{T,4}   # (Nr, Nr, Nμ, Nμ)
end

@inline HistR12R23Mu12Mu13(Nr::Integer, Nμ::Integer; T::Type{<:Integer}=Int) =
    HistR12R23Mu12Mu13(zeros(T, Nr, Nr, Nμ, Nμ))


# =========================== Binning =============================

@inline function bin_triangle_pre!(Hh::AbstractArray{<:Integer,4},
                                   r_a::Float64, r_b::Float64, μ_a::Float64, μ_b::Float64,
                                   rmin::Float64, invΔr::Float64, Nr::Int,
                                   invΔμ::Float64, Nμ::Int)
    ir = bin_index_inv(r_a, rmin, invΔr, Nr); ir == 0 && return
    jr = bin_index_inv(r_b, rmin, invΔr, Nr); jr == 0 && return
    im = bin_index_inv(μ_a, 0.0,  invΔμ, Nμ); im == 0 && return
    jm = bin_index_inv(μ_b, 0.0,  invΔμ, Nμ); jm == 0 && return
    @inbounds Hh[ir, jr, im, jm] += 1
end

# =========================== Two shortest ========================

@inline function two_shortest_using_r2(r2_12, r2_23, r2_13,
                                       md_12, md_23, md_13)
    a_r2, a_md = r2_12, md_12
    b_r2, b_md = r2_23, md_23
    c_r2, c_md = r2_13, md_13
    if b_r2 < a_r2
        a_r2, b_r2 = b_r2, a_r2
        a_md, b_md = b_md, a_md
    end
    if c_r2 < a_r2
        a_r2, c_r2 = c_r2, a_r2
        a_md, c_md = c_md, a_md
    end
    if c_r2 < b_r2
        return sqrt(a_r2), sqrt(c_r2), a_md, c_md
    else
        return sqrt(a_r2), sqrt(b_r2), a_md, b_md
    end
end

# =========================== Grid types ==========================

struct GridNonPeriodic
    x0::Float64; y0::Float64; z0::Float64
    cs::Float64              # cell size
    invcs::Float64
    nx::Int; ny::Int; nz::Int
    # CSR-like layout for points per cell:
    # starts[c]..(starts[c]+counts[c]-1) are indices into "pts"
    starts::Vector{Int}
    counts::Vector{Int}
    pts::Vector{Int}
end

struct GridPeriodic
    Lx::Float64; Ly::Float64; Lz::Float64
    hx::Float64; hy::Float64; hz::Float64 # half-lengths
    cs::Float64
    invcs::Float64
    nx::Int; ny::Int; nz::Int
    starts::Vector{Int}
    counts::Vector{Int}
    pts::Vector{Int}
end

# =========================== Build grids =========================

# Non-periodic: grid bounds from data min/max (padded slightly)
function build_grid_nonperiodic(X::AbstractVector{<:Real},
                                Y::AbstractVector{<:Real},
                                Z::AbstractVector{<:Real},
                                cellsize::Real)
    N = length(X)
    cs = Float64(cellsize)
    invcs = 1.0 / cs

    xmin = minimum(Float64.(X)); xmax = maximum(Float64.(X))
    ymin = minimum(Float64.(Y)); ymax = maximum(Float64.(Y))
    zmin = minimum(Float64.(Z)); zmax = maximum(Float64.(Z))
    # pad by tiny epsilon to avoid edge issues
    pad = 1e-9 * max(xmax - xmin, max(ymax - ymin, zmax - zmin))
    x0 = xmin - pad; y0 = ymin - pad; z0 = zmin - pad

    nx = max(1, Int(ceil((xmax - x0 + pad) * invcs)))
    ny = max(1, Int(ceil((ymax - y0 + pad) * invcs)))
    nz = max(1, Int(ceil((zmax - z0 + pad) * invcs)))
    ncells = nx * ny * nz

    counts = zeros(Int, ncells)
    # first pass: counts
    @inbounds for i in 1:N
        ix = Int(floor((Float64(X[i]) - x0) * invcs)); (ix < 0 || ix >= nx) && continue
        iy = Int(floor((Float64(Y[i]) - y0) * invcs)); (iy < 0 || iy >= ny) && continue
        iz = Int(floor((Float64(Z[i]) - z0) * invcs)); (iz < 0 || iz >= nz) && continue
        cid = ((ix * ny) + iy) * nz + iz + 1
        counts[cid] += 1
    end
    starts = similar(counts)
    s = 1
    @inbounds for c in 1:ncells
        starts[c] = s
        s += counts[c]
    end
    pts = Vector{Int}(undef, N)
    # second pass: fill
    fill!(counts, 0)  # reuse as write cursors
    @inbounds for i in 1:N
        ix = Int(floor((Float64(X[i]) - x0) * invcs)); (ix < 0 || ix >= nx) && continue
        iy = Int(floor((Float64(Y[i]) - y0) * invcs)); (iy < 0 || iy >= ny) && continue
        iz = Int(floor((Float64(Z[i]) - z0) * invcs)); (iz < 0 || iz >= nz) && continue
        cid = ((ix * ny) + iy) * nz + iz + 1
        pos = starts[cid] + counts[cid]
        pts[pos] = i
        counts[cid] += 1
    end
    return GridNonPeriodic(x0, y0, z0, cs, invcs, nx, ny, nz, starts, counts, pts)
end

# Periodic: cells tile the box exactly
function build_grid_periodic(X::AbstractVector{<:Real},
                             Y::AbstractVector{<:Real},
                             Z::AbstractVector{<:Real};
                             Lx::Real, Ly::Real, Lz::Real,
                             cellsize::Real)
    N = length(X)
    Lxf = Float64(Lx); Lyf = Float64(Ly); Lzf = Float64(Lz)
    cs_guess = Float64(cellsize)
    # choose integer cells per axis, at least 1
    nx = max(1, Int(floor(Lxf / cs_guess)))
    ny = max(1, Int(floor(Lyf / cs_guess)))
    nz = max(1, Int(floor(Lzf / cs_guess)))
    # recompute exact cell size so nx*cs = L
    csx = Lxf / nx; csy = Lyf / ny; csz = Lzf / nz
    # we’ll use axis-specific cell sizes; store the geometric mean as cs, but keep inv per axis locally
    # For mapping cells we’ll use separate invcs per axis inside counting routine.
    cs = (csx + csy + csz) / 3
    invcsx = 1.0 / csx; invcsy = 1.0 / csy; invcsz = 1.0 / csz
    ncells = nx * ny * nz

    counts = zeros(Int, ncells)
    # map to cell with periodic wrap
    @inline cellwrap(u, L, invcs, n) = begin
        k = Int(floor(u * invcs)) % n
        (k < 0) && (k += n)
        k
    end

    # counts
    @inbounds for i in 1:N
        ix = cellwrap(Float64(X[i]), Lxf, invcsx, nx)
        iy = cellwrap(Float64(Y[i]), Lyf, invcsy, ny)
        iz = cellwrap(Float64(Z[i]), Lzf, invcsz, nz)
        cid = ((ix * ny) + iy) * nz + iz + 1
        counts[cid] += 1
    end
    starts = similar(counts)
    s = 1
    @inbounds for c in 1:ncells
        starts[c] = s
        s += counts[c]
    end
    pts = Vector{Int}(undef, N)
    fill!(counts, 0)
    @inbounds for i in 1:N
        ix = cellwrap(Float64(X[i]), Lxf, invcsx, nx)
        iy = cellwrap(Float64(Y[i]), Lyf, invcsy, ny)
        iz = cellwrap(Float64(Z[i]), Lzf, invcsz, nz)
        cid = ((ix * ny) + iy) * nz + iz + 1
        pos = starts[cid] + counts[cid]
        pts[pos] = i
        counts[cid] += 1
    end
    return GridPeriodic(Lxf, Lyf, Lzf, 0.5Lxf, 0.5Lyf, 0.5Lzf, cs, 1.0/cs, nx, ny, nz, starts, counts, pts)
end

# ===================== Neighbor iteration ========================

# iterate all points in a cell (by linear id)
@inline function cell_points(starts::Vector{Int}, counts::Vector{Int}, pts::Vector{Int}, cid::Int)
    s = starts[cid]; c = counts[cid]
    return @view pts[s : s + c - 1]
end

# ======================= Counting kernels ========================

"""
Non-periodic, true LOS, grid-accelerated.

Inputs:
  X, Y, Z::AbstractVector{<:Real}
  rmin, rmax::Real, Nr::Integer, μmax::Real, Nμ::Integer
  cellsize::Real   (suggest: ≈ rmax)
"""
function count_triangles_grid!(X::AbstractVector{<:Real},
                               Y::AbstractVector{<:Real},
                               Z::AbstractVector{<:Real};
                               rmin::Real, rmax::Real, Nr::Integer,
                               μmax::Real, Nμ::Integer=2,
                               cellsize::Real=rmax)

    N = length(X)
    @assert length(Y) == N && length(Z) == N

    # constants
    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf * μmaxf
    NrI = Int(Nr); NμI = Int(Nμ)

    invΔr = NrI / (rmaxf - rminf)
    invΔμ = NμI / μmaxf
    r2min = rminf * rminf
    r2max = rmaxf * rmaxf

    # build grid
    G = build_grid_nonperiodic(X, Y, Z, cellsize)

    # per-thread histograms
    Hs = [zeros(Int, NrI, NrI, NμI, NμI) for _ in 1:Threads.nthreads()]

    # Thread-local neighbor buffer to avoid allocations
    buffers = [Vector{Int}(undef, 0) for _ in 1:Threads.nthreads()]
    foreach(b -> sizehint!(b, 2048), buffers)

    @inbounds Threads.@threads for ci in 1:(G.nx*G.ny*G.nz)
        tid = Threads.threadid()
        Hh  = Hs[tid]
        buf = buffers[tid]
        empty!(buf)

        # Decode cell coordinates
        ix = (ci - 1) ÷ (G.ny * G.nz)
        rem1 = (ci - 1) % (G.ny * G.nz)
        iy = rem1 ÷ G.nz
        iz = rem1 % G.nz

        pts_i = cell_points(G.starts, G.counts, G.pts, ci)
        length(pts_i) == 0 && continue

        # neighbor cell range (clamped)
        x0 = max(ix - 1, 0); x1 = min(ix + 1, G.nx - 1)
        y0 = max(iy - 1, 0); y1 = min(iy + 1, G.ny - 1)
        z0 = max(iz - 1, 0); z1 = min(iz + 1, G.nz - 1)

        for idx_i in pts_i
            xi = Float64(X[idx_i]); yi = Float64(Y[idx_i]); zi = Float64(Z[idx_i])

            # Build neighbor list j>i that pass (rmin..rmax) and μ threshold vs i
            empty!(buf)
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
                    push!(buf, j)
                end
            end

            # Now form triangles from neighbor list (j < k; both > i)
            m = length(buf)
            if m >= 2
                for a in 1:m-1
                    j = buf[a]
                    xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
                    dxij = xj - xi; dyij = yj - yi; dzij = zj - zi
                    r2_ij = dxij*dxij + dyij*dyij + dzij*dzij
                    md_ij = mu_true_los_numden2(xi, yi, zi, xj, yj, zj)  # reuse

                    for b in a+1:m
                        k = buf[b]
                        xk = Float64(X[k]); yk = Float64(Y[k]); zk = Float64(Z[k])

                        # i-k was pruned already (since k in buf), so only check j-k
                        dxjk = xk - xj; dyjk = yk - yj; dzjk = zk - zj
                        r2_jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk
                        (r2_jk < r2min || r2_jk >= r2max) && continue
                        md_jk = mu_true_los_numden2(xj, yj, zj, xk, yk, zk)
                        num2, den2 = md_jk
                        !(den2 == 0.0 || num2 < μmax2 * den2) && continue

                        dxik = xi - xk; dyik = yi - yk; dzik = zi - zk
                        r2_ik = dxik*dxik + dyik*dyik + dzik*dzik
                        # (i-k already passed earlier, so no μ check here)

                        r_a, r_b, md_a, md_b = two_shortest_using_r2(r2_ij, r2_jk, r2_ik,
                                                                     md_ij, md_jk,
                                                                     mu_true_los_numden2(xi, yi, zi, xk, yk, zk))
                        μ_a = mu_from_numden2(md_a...)
                        μ_b = mu_from_numden2(md_b...)
                        bin_triangle_pre!(Hh, r_a, r_b, μ_a, μ_b, rminf, invΔr, NrI, invΔμ, NμI)
                    end
                end
            end
        end
    end

    # reduce histograms
    H = HistR12R23Mu12Mu13(NrI, NμI; T=Int)
    @inbounds for t in 1:length(Hs)
        H.h .+= Hs[t]
    end
    return H
end

"""
Periodic, min-image separations, fixed z-LOS, grid-accelerated.
The box is [0,Lx)×[0,Ly)×[0,Lz).
"""
function count_triangles_periodic_grid!(X::AbstractVector{<:Real},
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

    # build periodic grid
    G = build_grid_periodic(X, Y, Z; Lx=Lx, Ly=Ly, Lz=Lz, cellsize=cellsize)
    Lxf, Lyf, Lzf = G.Lx, G.Ly, G.Lz
    hx, hy, hz = G.hx, G.hy, G.hz

    # axis-specific cell sizes for wrap math
    csx = Lxf / G.nx; csy = Lyf / G.ny; csz = Lzf / G.nz

    # per-thread histograms & buffers
    Hs = [zeros(Int, NrI, NrI, NμI, NμI) for _ in 1:Threads.nthreads()]
    buffers = [Vector{Int}(undef, 0) for _ in 1:Threads.nthreads()]
    foreach(b -> sizehint!(b, 2048), buffers)

    # branch-only min-image (fast)
    @inline function minimg_d(Δ::Float64, halfL::Float64, L::Float64)
        if Δ >  halfL
            Δ - L
        elseif Δ <= -halfL
            Δ + L
        else
            Δ
        end
    end

    # neighbor cell ranges (with periodic wrap)
    @inline function wrap(a::Int, n::Int)
        v = a % n
        (v < 0) && (v += n)
        return v
    end

    @inbounds Threads.@threads for ci in 1:(G.nx*G.ny*G.nz)
        tid = Threads.threadid()
        Hh  = Hs[tid]
        buf = buffers[tid]
        empty!(buf)

        ix = (ci - 1) ÷ (G.ny * G.nz)
        rem1 = (ci - 1) % (G.ny * G.nz)
        iy = rem1 ÷ G.nz
        iz = rem1 % G.nz

        pts_i = cell_points(G.starts, G.counts, G.pts, ci)
        length(pts_i) == 0 && continue

        for idx_i in pts_i
            xi = Float64(X[idx_i]); yi = Float64(Y[idx_i]); zi = Float64(Z[idx_i])

            # Build neighbor list of j>i around (ix,iy,iz) with periodic wrap
            empty!(buf)
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
                    push!(buf, j)
                end
            end

            # triangles from neighbor list
            m = length(buf)
            if m >= 2
                for a in 1:m-1
                    j = buf[a]
                    xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
                    dxij = minimg_d(xj - xi, hx, Lxf)
                    dyij = minimg_d(yj - yi, hy, Lyf)
                    dzij = minimg_d(zj - zi, hz, Lzf)
                    r2_ij = dxij*dxij + dyij*dyij + dzij*dzij
                    md_ij = mu_z_los_numden2(dxij, dyij, dzij)

                    for b in a+1:m
                        k = buf[b]
                        xk = Float64(X[k]); yk = Float64(Y[k]); zk = Float64(Z[k])

                        dxjk = minimg_d(xk - xj, hx, Lxf)
                        dyjk = minimg_d(yk - yj, hy, Lyf)
                        dzjk = minimg_d(zk - zj, hz, Lzf)
                        r2_jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk
                        (r2_jk < r2min || r2_jk >= r2max) && continue
                        md_jk = mu_z_los_numden2(dxjk, dyjk, dzjk)
                        num2, den2 = md_jk
                        !(den2 == 0.0 || num2 < μmax2 * den2) && continue

                        dxik = minimg_d(xi - xk, hx, Lxf)
                        dyik = minimg_d(yi - yk, hy, Lyf)
                        dzik = minimg_d(zi - zk, hz, Lzf)
                        r2_ik = dxik*dxik + dyik*dyik + dzik*dzik
                        # i-k μ already passed when k was added to buf

                        r_a, r_b, md_a, md_b =
                            two_shortest_using_r2(r2_ij, r2_jk, r2_ik, md_ij, md_jk, mu_z_los_numden2(dxik, dyik, dzik))
                        μ_a = mu_from_numden2(md_a...)
                        μ_b = mu_from_numden2(md_b...)
                        bin_triangle_pre!(Hh, r_a, r_b, μ_a, μ_b, rminf, invΔr, NrI, invΔμ, NμI)
                    end
                end
            end
        end
    end

    H = HistR12R23Mu12Mu13(NrI, NμI; T=Int)
    @inbounds for t in 1:length(Hs)
        H.h .+= Hs[t]
    end
    return H
end

