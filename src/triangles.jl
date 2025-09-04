using Base.Threads

# 4D joint histogram over (r12, r23, μ12, μ13)
mutable struct HistR12R23Mu12Mu13{T<:Integer}
    h::Array{T,4}   # size (Nr, Nr, Nμ, Nμ)
end
@inline HistR12R23Mu12Mu13(Nr::Integer, Nμ::Integer; T::Type{<:Integer}=Int) =
    HistR12R23Mu12Mu13(zeros(T, Nr, Nr, Nμ, Nμ))

# --- helper: deterministic 2-way swap used in tiny 3-item sort ---
@inline function cswap!(r_a::Int, r_b::Int, m_a::Int, m_b::Int, id_a::Int, id_b::Int)
    # Order by r-bin first, then μ-bin, then a fixed pair-id (1:(12), 2:(23), 3:(13))
    if (r_a > r_b) || (r_a == r_b && (m_a > m_b || (m_a == m_b && id_a > id_b)))
        r_a, r_b = r_b, r_a
        m_a, m_b = m_b, m_a
        id_a, id_b = id_b, id_a
    end
    return r_a, r_b, m_a, m_b, id_a, id_b
end


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

    b = TripletBinner(rmin, rmax, Nr, μmax, Nμ)
    g = build_grid(X, Y, Z, cellsize)

    nt = nthreads()
    Ht = [zeros(Int, b.Nr, b.Nr, b.Nμ, b.Nμ) for _ in 1:nt]

    # per-thread neighbor scratch for (i,j) that passed cuts (+ cached bins)
    Jidx = [Int[] for _ in 1:nt]
    Jr   = [Int[] for _ in 1:nt]   # r12 bin
    Jm   = [Int[] for _ in 1:nt]   # μ12 bin
    @inbounds for t in 1:nt
        sizehint!(Jidx[t], 256); sizehint!(Jr[t], 256); sizehint!(Jm[t], 256)
    end

    @threads for i in 1:N
        tid = threadid()
        H    = Ht[tid]
        jidx = Jidx[tid]; jr12 = Jr[tid]; jm12 = Jm[tid]
        empty!(jidx); empty!(jr12); empty!(jm12)

        x1 = float(X[i]); y1 = float(Y[i]); z1 = float(Z[i])
        cxi, cyi, czi = cell_of(g, x1, y1, z1)

        # 1) gather neighbors (i,j)
        @inbounds for dz in -1:1, dy in -1:1, dx in -1:1
            cx = cxi + dx; cy = cyi + dy; cz = czi + dz
            if (0 <= cx < g.nx) & (0 <= cy < g.ny) & (0 <= cz < g.nz)
                cid = cell_id(g, cx, cy, cz)
                s = g.offs[cid]; e = g.offs[cid+1]-1
                for p in s:e
                    j = g.idxs[p]
                    j <= i && continue

                    x2 = float(X[j]); y2 = float(Y[j]); z2 = float(Z[j])

                    r12_sq = dist2(x1,y1,z1, x2,y2,z2)
                    (r12_sq <= b.rmin2 || r12_sq >= b.rmax2) && continue

                    μ12_sq = mu2_true_los(x1,y1,z1, x2,y2,z2, r12_sq)
                    μ12_sq >= b.μmax2 && continue

                    r12b = Int(floor((sqrt(r12_sq) - b.rmin) * b.invΔr)) + 1
                    μ12b = Int(floor( sqrt(μ12_sq)            * b.invΔμ)) + 1
                    (r12b < 1 || r12b > b.Nr || μ12b < 1 || μ12b > b.Nμ) && continue

                    push!(jidx, j); push!(jr12, r12b); push!(jm12, μ12b)
                end
            end
        end

        # 2) form (j,k), enforce (j,k) cuts, then sort (r,μ) and bin (short, medium)
        L = length(jidx)
        @inbounds for a in 1:(L-1)
            j = jidx[a]
            x2 = float(X[j]); y2 = float(Y[j]); z2 = float(Z[j])
            r12b = jr12[a]; μ12b = jm12[a]

            for b2 in (a+1):L
                k = jidx[b2]

                r23_sq = dist2(x2,y2,z2, float(X[k]), float(Y[k]), float(Z[k]))
                (r23_sq <= b.rmin2 || r23_sq >= b.rmax2) && continue
                μ23_sq = mu2_true_los(x2,y2,z2, float(X[k]),float(Y[k]),float(Z[k]), r23_sq)
                μ23_sq >= b.μmax2 && continue

                r23b = Int(floor((sqrt(r23_sq) - b.rmin) * b.invΔr)) + 1
                (r23b < 1 || r23b > b.Nr) && continue
                μ23b = Int(floor(sqrt(μ23_sq) * b.invΔμ)) + 1
                (μ23b < 1 || μ23b > b.Nμ) && continue

                # (i,k) from cache
                r13b = jr12[b2]; μ13b = jm12[b2]

                # sort by r-bin (short→long), carry μ; tie-break by μ then pair-id
                r1, r2, r3 = r12b, r23b, r13b
                m1, m2, m3 = μ12b, μ23b, μ13b
                id1, id2, id3 = 1, 2, 3

                r1, r2, m1, m2, id1, id2 = cswap!(r1, r2, m1, m2, id1, id2)
                r2, r3, m2, m3, id2, id3 = cswap!(r2, r3, m2, m3, id2, id3)
                r1, r2, m1, m2, id1, id2 = cswap!(r1, r2, m1, m2, id1, id2)

                H[r1, r2, m1, m2] += 1
            end
        end
    end

    H = HistR12R23Mu12Mu13(Nr, Nμ)
    @inbounds for t in 1:nt
        H.h .+= Ht[t]
    end
    return H
end


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

    Jidx = [Int[] for _ in 1:nt]
    Jr   = [Int[] for _ in 1:nt]
    Jm   = [Int[] for _ in 1:nt]
    @inbounds for t in 1:nt
        sizehint!(Jidx[t], 256); sizehint!(Jr[t], 256); sizehint!(Jm[t], 256)
    end

    @threads for i in 1:N
        tid = threadid()
        H    = Ht[tid]
        jidx = Jidx[tid]; jr12 = Jr[tid]; jm12 = Jm[tid]
        empty!(jidx); empty!(jr12); empty!(jm12)

        x1 = float(X[i]); y1 = float(Y[i]); z1 = float(Z[i])
        cxi, cyi, czi = cell_of(g, x1, y1, z1)

        # 1) periodic neighbors (min-image), fixed 3×3×3 stencil
        @inbounds for dz in -1:1, dy in -1:1, dx in -1:1
            cx = wrap(cxi + dx, g.nx); cy = wrap(cyi + dy, g.ny); cz = wrap(czi + dz, g.nz)
            cid = cell_id(g, cx, cy, cz)
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

        # 2) (j,k) + sorting → bin (short, medium)
        L = length(jidx)
        @inbounds for a in 1:(L-1)
            j = jidx[a]
            x2 = float(X[j]); y2 = float(Y[j]); z2 = float(Z[j])
            r12b = jr12[a]; μ12b = jm12[a]

            for b2 in (a+1):L
                k  = jidx[b2]
                xk = float(X[k]); yk = float(Y[k]); zk = float(Z[k])

                r23_sq, dx23, dy23, dz23 = dist2_periodic(x2,y2,z2, xk,yk,zk, Lx,Ly,Lz)
                (r23_sq <= b.rmin2 || r23_sq >= b.rmax2) && continue

                μ23_sq = mu2_zlos_from_deltas(dx23, dy23, dz23, r23_sq)
                μ23_sq >= b.μmax2 && continue

                r23b = Int(floor((sqrt(r23_sq) - b.rmin) * b.invΔr)) + 1
                (r23b < 1 || r23b > b.Nr) && continue
                μ23b = Int(floor(sqrt(μ23_sq) * b.invΔμ)) + 1
                (μ23b < 1 || μ23b > b.Nμ) && continue

                r13b = jr12[b2]; μ13b = jm12[b2]   # cached (i,k)

                r1, r2, r3 = r12b, r23b, r13b
                m1, m2, m3 = μ12b, μ23b, μ13b
                id1, id2, id3 = 1, 2, 3

                r1, r2, m1, m2, id1, id2 = cswap!(r1, r2, m1, m2, id1, id2)
                r2, r3, m2, m3, id2, id3 = cswap!(r2, r3, m2, m3, id2, id3)
                r1, r2, m1, m2, id1, id2 = cswap!(r1, r2, m1, m2, id1, id2)

                H[r1, r2, m1, m2] += 1
            end
        end
    end

    H = HistR12R23Mu12Mu13(Nr, Nμ)
    @inbounds for t in 1:nt
        H.h .+= Ht[t]
    end
    return H
end

