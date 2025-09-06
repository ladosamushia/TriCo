# src/triangles_mixed.jl

using Base.Threads

############ Catalog wrapper ############

"""
Lightweight catalog wrapper for mixed-triangle counting.
Holds references to X,Y,Z (no copies).
"""
struct TriCat{T<:AbstractVector}
    X::T; Y::T; Z::T
end
TriCat(X::AbstractVector, Y::AbstractVector, Z::AbstractVector) = TriCat{typeof(X)}(X,Y,Z)

@inline _same_catalog(A::TriCat, B::TriCat) =
    (A.X === B.X) & (A.Y === B.Y) & (A.Z === B.Z)

############ Public API ############

"""
Non-periodic mixed-catalog triangle counts with true LOS (pair midpoint).

Counts triangles with one vertex from each of A,B,C.
If two (or three) of A,B,C are the same catalog object, triangles are
deduplicated by enforcing index order only within the repeated catalogs:
- if A===B, require i<j
- if A===C, require i<k
- if B===C, require j<k

Sides are sorted by length before binning, so bins are label-free.
"""
function count_triangles_mixed!(A::TriCat, B::TriCat, C::TriCat;
    rmin::Real, rmax::Real, Nr::Integer,
    μmax::Real, Nμ::Integer=2, cellsize::Real=rmax)

    # constants (match triangles.jl style)
    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf * μmaxf
    NrI = Int(Nr); NμI = Int(Nμ)

    invΔr = NrI / (rmaxf - rminf)
    invΔμ = NμI / μmaxf
    r2min = rminf * rminf
    r2max = rmaxf * rmaxf

    # grids per catalog
    gA = build_grid_nonperiodic(A.X, A.Y, A.Z, cellsize)
    gB = build_grid_nonperiodic(B.X, B.Y, B.Z, cellsize)
    gC = build_grid_nonperiodic(C.X, C.Y, C.Z, cellsize)

    # per-thread histograms and one neighbor buffer for B
    Hs = [zeros(Int, NrI, NrI, NμI, NμI) for _ in 1:Threads.nthreads()]
    bufB = [Vector{Int}(undef, 0) for _ in 1:Threads.nthreads()]
    foreach(b -> sizehint!(b, 2048), bufB)

    sameAB = _same_catalog(A,B)
    sameAC = _same_catalog(A,C)
    sameBC = _same_catalog(B,C)

    # Map world coordinate → (ix,iy,iz) in a nonperiodic grid (clamped)
    @inline function grid_coords_nonperp(x::Float64, y::Float64, z::Float64, G::GridNonPeriodic)
        ix = Int(floor((x - G.x0) * G.invcs)); ix = ifelse(ix < 0, 0, ifelse(ix >= G.nx, G.nx-1, ix))
        iy = Int(floor((y - G.y0) * G.invcs)); iy = ifelse(iy < 0, 0, ifelse(iy >= G.ny, G.ny-1, iy))
        iz = Int(floor((z - G.z0) * G.invcs)); iz = ifelse(iz < 0, 0, ifelse(iz >= G.nz, G.nz-1, iz))
        return ix,iy,iz
    end

    @inbounds Threads.@threads for ci in 1:(gA.nx*gA.ny*gA.nz)
        tid = Threads.threadid()
        Hh  = Hs[tid]
        bB  = bufB[tid]
        empty!(bB)

        # A-cell coords
        aix = (ci - 1) ÷ (gA.ny * gA.nz)
        rem1 = (ci - 1) % (gA.ny * gA.nz)
        aiy = rem1 ÷ gA.nz
        aiz = rem1 % gA.nz

        ptsA = cell_points(gA.starts, gA.counts, gA.pts, ci)
        length(ptsA) == 0 && continue

        # A neighbor cell range (clamped)
        ax0 = max(aix - 1, 0); ax1 = min(aix + 1, gA.nx - 1)
        ay0 = max(aiy - 1, 0); ay1 = min(aiy + 1, gA.ny - 1)
        az0 = max(aiz - 1, 0); az1 = min(aiz + 1, gA.nz - 1)

        for i in ptsA
            xi = Float64(A.X[i]); yi = Float64(A.Y[i]); zi = Float64(A.Z[i])

            # Determine B-cell neighborhood around the *position of i* in B's grid
            bix, biy, biz = grid_coords_nonperp(xi, yi, zi, gB)
            bx0 = max(bix - 1, 0); bx1 = min(bix + 1, gB.nx - 1)
            by0 = max(biy - 1, 0); by1 = min(biy + 1, gB.ny - 1)
            bz0 = max(biz - 1, 0); bz1 = min(biz + 1, gB.nz - 1)

            # Build neighbor list of j (from B) that pass r and μ with i
            empty!(bB)
            @inbounds for bx in bx0:bx1, by in by0:by1, bz in bz0:bz1
                bcid = ((bx * gB.ny) + by) * gB.nz + bz + 1
                for j in cell_points(gB.starts, gB.counts, gB.pts, bcid)
                    if sameAB && !(i < j); continue; end  # dedup within repeated catalog
                    xj = Float64(B.X[j]); yj = Float64(B.Y[j]); zj = Float64(B.Z[j])
                    dx = xj - xi; dy = yj - yi; dz = zj - zi
                    r2 = dx*dx + dy*dy + dz*dz
                    (r2 < r2min || r2 >= r2max) && continue
                    num2, den2 = mu_true_los_numden2(xi, yi, zi, xj, yj, zj)
                    !(den2 == 0.0 || num2 < μmax2 * den2) && continue
                    push!(bB, j)
                end
            end

            m = length(bB)
            m < 1 && continue

            # For each accepted (i∈A, j∈B), try k∈C near i
            # Neighborhood in C around i:
            cix, ciy, ciz = grid_coords_nonperp(xi, yi, zi, gC)
            cx0 = max(cix - 1, 0); cx1 = min(cix + 1, gC.nx - 1)
            cy0 = max(ciy - 1, 0); cy1 = min(ciy + 1, gC.ny - 1)
            cz0 = max(ciz - 1, 0); cz1 = min(ciz + 1, gC.nz - 1)

            @inbounds for a in 1:m
                j = bB[a]
                xj = Float64(B.X[j]); yj = Float64(B.Y[j]); zj = Float64(B.Z[j])
                dxij = xj - xi; dyij = yj - yi; dzij = zj - zi
                r2_ij = dxij*dxij + dyij*dyij + dzij*dzij
                md_ij = mu_true_los_numden2(xi, yi, zi, xj, yj, zj)

                for cx in cx0:cx1, cy in cy0:cy1, cz in cz0:cz1
                    ccid = ((cx * gC.ny) + cy) * gC.nz + cz + 1
                    for k in cell_points(gC.starts, gC.counts, gC.pts, ccid)
                        if sameAC && !(i < k); continue; end
                        if sameBC && !(j < k); continue; end

                        xk = Float64(C.X[k]); yk = Float64(C.Y[k]); zk = Float64(C.Z[k])

                        # distances
                        dxik = xk - xi; dyik = yk - yi; dzik = zk - zi
                        rik2 = dxik*dxik + dyik*dyik + dzik*dzik
                        (rik2 < r2min || rik2 >= r2max) && continue

                        dxjk = xk - xj; dyjk = yk - yj; dzjk = zk - zj
                        rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk
                        (rjk2 < r2min || rjk2 >= r2max) && continue

                        # μ checks
                        md_ik = mu_true_los_numden2(xi, yi, zi, xk, yk, zk)
                        num2, den2 = md_ik
                        !(den2 == 0.0 || num2 < μmax2 * den2) && continue
                        md_jk = mu_true_los_numden2(xj, yj, zj, xk, yk, zk)
                        num2, den2 = md_jk
                        !(den2 == 0.0 || num2 < μmax2 * den2) && continue

                        # Canonicalize by side lengths (label-free) and bin
                        r_a, r_b, md_a, md_b =
                            two_shortest_using_r2(r2_ij, rjk2, rik2, md_ij, md_jk, md_ik)
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
Periodic mixed-catalog triangles (min-image, z-LOS).
Box is [0,Lx)×[0,Ly)×[0,Lz).
"""
function count_triangles_mixed_periodic!(A::TriCat, B::TriCat, C::TriCat;
    Lx::Real, Ly::Real, Lz::Real,
    rmin::Real, rmax::Real, Nr::Integer,
    μmax::Real, Nμ::Integer=2, cellsize::Real=rmax)

    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf * μmaxf
    NrI = Int(Nr); NμI = Int(Nμ)

    invΔr = NrI / (rmaxf - rminf)
    invΔμ = NμI / μmaxf
    r2min = rminf * rminf
    r2max = rmaxf * rmaxf

    gA = build_grid_periodic(A.X, A.Y, A.Z; Lx=Lx, Ly=Ly, Lz=Lz, cellsize=cellsize)
    gB = build_grid_periodic(B.X, B.Y, B.Z; Lx=Lx, Ly=Ly, Lz=Lz, cellsize=cellsize)
    gC = build_grid_periodic(C.X, C.Y, C.Z; Lx=Lx, Ly=Ly, Lz=Lz, cellsize=cellsize)

    Lxf, Lyf, Lzf = gA.Lx, gA.Ly, gA.Lz
    hx, hy, hz = gA.hx, gA.hy, gA.hz

    # per-thread hist + B-neighbors buffer
    Hs = [zeros(Int, NrI, NrI, NμI, NμI) for _ in 1:Threads.nthreads()]
    bufB = [Vector{Int}(undef, 0) for _ in 1:Threads.nthreads()]
    foreach(b -> sizehint!(b, 2048), bufB)

    sameAB = _same_catalog(A,B)
    sameAC = _same_catalog(A,C)
    sameBC = _same_catalog(B,C)

    # min-image and grid wrapping helpers (mirror triangles.jl style)
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

    # Map world coordinate → (ix,iy,iz) in a periodic grid (wrapped)
    @inline function grid_coords_perp(x::Float64, y::Float64, z::Float64, G::GridPeriodic)
        csx = G.Lx / G.nx; csy = G.Ly / G.ny; csz = G.Lz / G.nz
        ix = Int(floor(x / csx)) % G.nx; ix < 0 && (ix += G.nx)
        iy = Int(floor(y / csy)) % G.ny; iy < 0 && (iy += G.ny)
        iz = Int(floor(z / csz)) % G.nz; iz < 0 && (iz += G.nz)
        return ix,iy,iz
    end

    @inbounds Threads.@threads for ci in 1:(gA.nx*gA.ny*gA.nz)
        tid = Threads.threadid()
        Hh  = Hs[tid]
        bB  = bufB[tid]
        empty!(bB)

        aix = (ci - 1) ÷ (gA.ny * gA.nz)
        rem1 = (ci - 1) % (gA.ny * gA.nz)
        aiy = rem1 ÷ gA.nz
        aiz = rem1 % gA.nz

        ptsA = cell_points(gA.starts, gA.counts, gA.pts, ci)
        length(ptsA) == 0 && continue

        for i in ptsA
            xi = Float64(A.X[i]); yi = Float64(A.Y[i]); zi = Float64(A.Z[i])

            # B neighbors around i's position
            bix, biy, biz = grid_coords_perp(xi, yi, zi, gB)
            @inbounds begin
                empty!(bB)
                for dx in -1:1, dy in -1:1, dz in -1:1
                    jx = wrap(bix + dx, gB.nx)
                    jy = wrap(biy + dy, gB.ny)
                    jz = wrap(biz + dz, gB.nz)
                    bcid = ((jx * gB.ny) + jy) * gB.nz + jz + 1
                    for j in cell_points(gB.starts, gB.counts, gB.pts, bcid)
                        if sameAB && !(i < j); continue; end
                        xj = Float64(B.X[j]); yj = Float64(B.Y[j]); zj = Float64(B.Z[j])
                        dxij = minimg_d(xj - xi, hx, Lxf)
                        dyij = minimg_d(yj - yi, hy, Lyf)
                        dzij = minimg_d(zj - zi, hz, Lzf)
                        r2 = dxij*dxij + dyij*dyij + dzij*dzij
                        (r2 < r2min || r2 >= r2max) && continue
                        num2, den2 = mu_z_los_numden2(dxij, dyij, dzij)
                        !(den2 == 0.0 || num2 < μmax2 * den2) && continue
                        push!(bB, j)
                    end
                end
            end

            m = length(bB)
            m < 1 && continue

            # C neighbors around i's position
            cix, ciy, ciz = grid_coords_perp(xi, yi, zi, gC)

            @inbounds for a in 1:m
                j = bB[a]
                xj = Float64(B.X[j]); yj = Float64(B.Y[j]); zj = Float64(B.Z[j])
                dxij = minimg_d(xj - xi, hx, Lxf)
                dyij = minimg_d(yj - yi, hy, Lyf)
                dzij = minimg_d(zj - zi, hz, Lzf)
                r2_ij = dxij*dxij + dyij*dyij + dzij*dzij
                md_ij = mu_z_los_numden2(dxij, dyij, dzij)

                for dx in -1:1, dy in -1:1, dz in -1:1
                    kx = wrap(cix + dx, gC.nx)
                    ky = wrap(ciy + dy, gC.ny)
                    kz = wrap(ciz + dz, gC.nz)
                    ccid = ((kx * gC.ny) + ky) * gC.nz + kz + 1
                    for k in cell_points(gC.starts, gC.counts, gC.pts, ccid)
                        if sameAC && !(i < k); continue; end
                        if sameBC && !(j < k); continue; end

                        xk = Float64(C.X[k]); yk = Float64(C.Y[k]); zk = Float64(C.Z[k])

                        dxik = minimg_d(xk - xi, hx, Lxf)
                        dyik = minimg_d(yk - yi, hy, Lyf)
                        dzik = minimg_d(zk - zi, hz, Lzf)
                        rik2 = dxik*dxik + dyik*dyik + dzik*dzik
                        (rik2 < r2min || rik2 >= r2max) && continue

                        dxjk = minimg_d(xk - xj, hx, Lxf)
                        dyjk = minimg_d(yk - yj, hy, Lyf)
                        dzjk = minimg_d(zk - zj, hz, Lzf)
                        rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk
                        (rjk2 < r2min || rjk2 >= r2max) && continue

                        md_ik = mu_z_los_numden2(dxik, dyik, dzik)
                        num2, den2 = md_ik
                        !(den2 == 0.0 || num2 < μmax2 * den2) && continue
                        md_jk = mu_z_los_numden2(dxjk, dyjk, dzjk)
                        num2, den2 = md_jk
                        !(den2 == 0.0 || num2 < μmax2 * den2) && continue

                        r_a, r_b, md_a, md_b =
                            two_shortest_using_r2(r2_ij, rjk2, rik2, md_ij, md_jk, md_ik)
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

