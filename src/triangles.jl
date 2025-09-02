using Base.Threads

# 4D joint histogram over (r12, r23, μ12, μ13)
mutable struct HistR12R23Mu12Mu13{T<:Integer}
    h::Array{T,4}   # size (Nr, Nr, Nμ, Nμ)
end
@inline HistR12R23Mu12Mu13(Nr::Integer, Nμ::Integer; T::Type{<:Integer}=Int) =
    HistR12R23Mu12Mu13(zeros(T, Nr, Nr, Nμ, Nμ))

"""
    count_triangles!(X, Y, Z; rmin, rmax, Nr, μmax, Nμ=2, cellsize=rmax)

Count triangles where **all pairs** satisfy: `rmin < r_ij < rmax` and `μ_ij < μmax`.
Return a **4D joint histogram** over `(r12, r23, μ12, μ13)` with shape `(Nr, Nr, Nμ, Nμ)`.

- Binning is **linear in r and μ** (half-open `[min, max)`).
- `μ23` is used **only for the cut**, not binned.
- Fast via uniform grid (cell size ≈ `rmax`) + per-thread histograms.

Assumes non-periodic Euclidean coords and clean inputs.
"""
function count_triangles!(X::AbstractVector{<:Real},
                          Y::AbstractVector{<:Real},
                          Z::AbstractVector{<:Real};
                          rmin::Real, rmax::Real, Nr::Integer,
                          μmax::Real, Nμ::Integer=2,
                          cellsize::Real=rmax)

    N = length(X); @assert length(Y) == N == length(Z)

    # Cached binning constants & fast range checks (defined in binning.jl)
    b = TripletBinner(rmin, rmax, Nr, μmax, Nμ)

    # Uniform grid for neighbor pruning (defined in grid.jl)
    g = build_grid(X, Y, Z, cellsize)

    # Per-thread 4D histograms to avoid locks
    nt = nthreads()
    Ht = [zeros(Int, b.Nr, b.Nr, b.Nμ, b.Nμ) for _ in 1:nt]

    # Per-thread scratch: neighbor list for i plus cached (r12-bin, μ12-bin)
    Jidx = [Int[] for _ in 1:nt]
    Jr   = [Int[] for _ in 1:nt]  # r12 bins
    Jm   = [Int[] for _ in 1:nt]  # μ12 bins
    @inbounds for t in 1:nt
        sizehint!(Jidx[t], 256); sizehint!(Jr[t], 256); sizehint!(Jm[t], 256)
    end

    @threads for i in 1:N
        tid = threadid()
        H = Ht[tid]
        jidx = Jidx[tid]; jr12 = Jr[tid]; jm12 = Jm[tid]
        empty!(jidx); empty!(jr12); empty!(jm12)

        x1 = float(X[i]); y1 = float(Y[i]); z1 = float(Z[i])
        cxi, cyi, czi = cell_of(g, x1, y1, z1)

        # 1) collect neighbors j of i that pass (r12, μ12) cuts; store (j, r12_bin, μ12_bin)
        @inbounds for dz in -1:1, dy in -1:1, dx in -1:1
            cx = cxi + dx; cy = cyi + dy; cz = czi + dz
            if (0 <= cx < g.nx) & (0 <= cy < g.ny) & (0 <= cz < g.nz)
                cid = cell_id(g, cx, cy, cz)
                s = g.offs[cid]; e = g.offs[cid+1] - 1
                for p in s:e
                    j = g.idxs[p]
                    j <= i && continue  # enforce i<j to avoid duplicates

                    x2 = float(X[j]); y2 = float(Y[j]); z2 = float(Z[j])

                    r12_sq = dist2(x1,y1,z1, x2,y2,z2)
                    (r12_sq <= b.rmin2 || r12_sq >= b.rmax2) && continue

                    μ12_sq = mu2_true_los(x1,y1,z1, x2,y2,z2, r12_sq)
                    μ12_sq >= b.μmax2 && continue

                    # uniform bins in linear r and μ
                    r12b = Int(floor((sqrt(r12_sq) - b.rmin) * b.invΔr)) + 1
                    μ12b = Int(floor( sqrt(μ12_sq)            * b.invΔμ)) + 1
                    (r12b < 1 || r12b > b.Nr || μ12b < 1 || μ12b > b.Nμ) && continue

                    push!(jidx, j); push!(jr12, r12b); push!(jm12, μ12b)
                end
            end
        end

        # 2) for each distinct neighbor pair (j,k), apply (j,k) cuts and bin r23 + cached μ13
        L = length(jidx)
        @inbounds for a in 1:(L-1)
            j = jidx[a]
            x2 = float(X[j]); y2 = float(Y[j]); z2 = float(Z[j])
            r12b = jr12[a]; μ12b = jm12[a]
            for b2 in (a+1):L
                k = jidx[b2]

                # (j,k) constraints: r23 and μ23
                r23_sq = dist2(x2,y2,z2, float(X[k]), float(Y[k]), float(Z[k]))
                (r23_sq <= b.rmin2 || r23_sq >= b.rmax2) && continue
                μ23_sq = mu2_true_los(x2,y2,z2, float(X[k]),float(Y[k]),float(Z[k]), r23_sq)
                μ23_sq >= b.μmax2 && continue

                # Passed all cuts → bin r23; μ13 bin is cached as jm12[b2]
                r23b = Int(floor((sqrt(r23_sq) - b.rmin) * b.invΔr)) + 1
                (r23b < 1 || r23b > b.Nr) && continue
                μ13b = jm12[b2]

                H[r12b, r23b, μ12b, μ13b] += 1
            end
        end
    end

    # Reduce per-thread histograms
    H = HistR12R23Mu12Mu13(Nr, Nμ)
    @inbounds for t in 1:nt
        H.h .+= Ht[t]
    end
    return H
end

using Base.Threads

"""
    count_triangles_periodic!(X,Y,Z; Lx,Ly,Lz, rmin,rmax,Nr, μmax, Nμ=2, cellsize=rmax)

Periodic-box version with **z as LOS** (μ² = Δz² / r²).

- Box lengths can differ by axis: `Lx`, `Ly`, `Lz`.
- Points are assumed to be in `[0,Lx)×[0,Ly)×[0,Lz)`.
- Cuts for **all three pairs**: `rmin < r_ij < rmax` and `μ_ij < μmax`.
- Returns a 4D joint histogram over `(r12, r23, μ12, μ13)` with size `(Nr, Nr, Nμ, Nμ)`.
"""
function count_triangles_periodic!(X::AbstractVector{<:Real},
                                        Y::AbstractVector{<:Real},
                                        Z::AbstractVector{<:Real};
                                        Lx::Real, Ly::Real, Lz::Real,
                                        rmin::Real, rmax::Real, Nr::Integer,
                                        μmax::Real, Nμ::Integer=2,
                                        cellsize::Real=rmax)

    N = length(X); @assert length(Y)==N==length(Z)
    Lx = float(Lx); Ly = float(Ly); Lz = float(Lz)

    b = TripletBinner(rmin, rmax, Nr, μmax, Nμ)
    g = build_grid_periodic(X, Y, Z, Lx, Ly, Lz, cellsize)

    nt = nthreads()
    Ht = [zeros(Int, b.Nr, b.Nr, b.Nμ, b.Nμ) for _ in 1:nt]

    # neighbor scratch: store j and their (r12, μ12) bins
    Jidx = [Int[] for _ in 1:nt]
    Jr   = [Int[] for _ in 1:nt]
    Jm   = [Int[] for _ in 1:nt]
    @inbounds for t in 1:nt
        sizehint!(Jidx[t], 256); sizehint!(Jr[t], 256); sizehint!(Jm[t], 256)
    end

    @threads for i in 1:N
        tid = threadid()
        H = Ht[tid]
        jidx = Jidx[tid]; jr12 = Jr[tid]; jm12 = Jm[tid]
        empty!(jidx); empty!(jr12); empty!(jm12)

        x1 = float(X[i]); y1 = float(Y[i]); z1 = float(Z[i])
        cxi, cyi, czi = cell_of(g, x1,y1,z1)

        # 1) Gather neighbors j of i that pass (r12, μ12) with periodic min-image distances
        @inbounds for dz in -1:1, dy in -1:1, dx in -1:1
            cx = wrap(cxi + dx, g.nx); cy = wrap(cyi + dy, g.ny); cz = wrap(czi + dz, g.nz)
            cid = cell_id(g, cx,cy,cz)
            s = g.offs[cid]; e = g.offs[cid+1]-1
            for p in s:e
                j = g.idxs[p]
                j <= i && continue
                x2 = float(X[j]); y2 = float(Y[j]); z2 = float(Z[j])

                r12_sq, dx12, dy12, dz12 = dist2_periodic(x1,y1,z1, x2,y2,z2, Lx,Ly,Lz)
                (r12_sq <= b.rmin2 || r12_sq >= b.rmax2) && continue

                μ12_sq = mu2_zlos_from_deltas(dx12, dy12, dz12, r12_sq)
                μ12_sq >= b.μmax2 && continue

                r12b = Int(floor((sqrt(r12_sq) - b.rmin) * b.invΔr)) + 1
                μ12b = Int(floor( sqrt(μ12_sq)            * b.invΔμ)) + 1
                (r12b < 1 || r12b > b.Nr || μ12b < 1 || μ12b > b.Nμ) && continue

                push!(jidx, j); push!(jr12, r12b); push!(jm12, μ12b)
            end
        end

        # 2) For each neighbor pair (j,k): enforce (j,k) cuts, then bin r23 and cached μ13
        L = length(jidx)
        @inbounds for a in 1:(L-1)
            j = jidx[a]
            x2 = float(X[j]); y2 = float(Y[j]); z2 = float(Z[j])
            r12b = jr12[a]; μ12b = jm12[a]
            for b2 in (a+1):L
                k = jidx[b2]

                r23_sq, dx23, dy23, dz23 = dist2_periodic(x2,y2,z2, float(X[k]),float(Y[k]),float(Z[k]), Lx,Ly,Lz)
                (r23_sq <= b.rmin2 || r23_sq >= b.rmax2) && continue
                μ23_sq = mu2_zlos_from_deltas(dx23, dy23, dz23, r23_sq)
                μ23_sq >= b.μmax2 && continue

                r23b = Int(floor((sqrt(r23_sq) - b.rmin) * b.invΔr)) + 1
                (r23b < 1 || r23b > b.Nr) && continue
                μ13b = jm12[b2]  # cached μ-bin for (i,k)

                H[r12b, r23b, μ12b, μ13b] += 1
            end
        end
    end

    H = HistR12R23Mu12Mu13(Nr, Nμ)
    @inbounds for t in 1:nt
        H.h .+= Ht[t]
    end
    return H
end

