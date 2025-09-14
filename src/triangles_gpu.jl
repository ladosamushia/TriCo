# src/triangles_gpu.jl
module TrianglesGPU

using CUDA
using StaticArrays
using ..PairsUtils: bin_index_inv
using ..: HistR12R23Mu12Mu13

# ------------------ Small GPU-safe helpers ------------------

# Return (num2, den2) for true LOS definition:
# μ^2 = ((s ⋅ (r_i + r_j))^2) / (|s|^2 * |r_i + r_j|^2)
@inline function gpu_mu_true_los_numden2(xi,yi,zi, xj,yj,zj)
    sx = xj - xi; sy = yj - yi; sz = zj - zi
    tx = (xi + xj); ty = (yi + yj); tz = (zi + zj)
    num = (sx*tx + sy*ty + sz*tz)
    num2 = num*num
    r2   = sx*sx + sy*sy + sz*sz
    t2   = tx*tx + ty*ty + tz*tz
    den2 = r2 * t2
    return num2, den2
end

# Periodic (fixed z LOS): μ^2 = Δz^2 / r^2 → (num2, den2) = (dz^2, r^2)
@inline function gpu_mu_z_los_numden2(dx,dy,dz)
    r2 = dx*dx + dy*dy + dz*dz
    return dz*dz, r2
end

# Compute μ from (num2, den2) safely (den2==0 ⇒ μ=0)
@inline function gpu_mu_from_numden2(num2, den2)
    den2 == 0.0 && return 0.0
    # Clamp tiny negatives from FP noise and sqrt the ratio
    v = num2 / den2
    v < 0.0 && (v = 0.0)
    return sqrt(v)
end

# Return the two shortest edges by r^2; with associated md tuples
@inline function gpu_two_shortest_using_r2(r2_12, r2_23, r2_13,
                                           md_12_num2, md_12_den2,
                                           md_23_num2, md_23_den2,
                                           md_13_num2, md_13_den2)
    a_r2 = r2_12; a_n=md_12_num2; a_d=md_12_den2
    b_r2 = r2_23; b_n=md_23_num2; b_d=md_23_den2
    c_r2 = r2_13; c_n=md_13_num2; c_d=md_13_den2
    if b_r2 < a_r2
        a_r2, b_r2 = b_r2, a_r2
        a_n,  b_n  = b_n,  a_n
        a_d,  b_d  = b_d,  a_d
    end
    if c_r2 < a_r2
        a_r2, c_r2 = c_r2, a_r2
        a_n,  c_n  = c_n,  a_n
        a_d,  c_d  = c_d,  a_d
    end
    if c_r2 < b_r2
        return sqrt(a_r2), sqrt(c_r2), a_n, a_d, c_n, c_d
    else
        return sqrt(a_r2), sqrt(b_r2), a_n, a_d, b_n, b_d
    end
end

# Device bin-indexer (inverse linear binning). Returns 0 if out-of-range.
@inline function gpu_bin_index_inv(x::Float64, xmin::Float64, invΔ::Float64, N::Int32)
    i = Int32(floor((x - xmin) * invΔ)) + 1
    return (i < 1 || i > N) ? Int32(0) : i
end

# ------------------ GPU grid mirrors ------------------

# Non-periodic grid mirror for GPU (CSR cell layout)
struct GPUGridNonPeriodic
    x0::Float64; y0::Float64; z0::Float64
    invcs::Float64
    nx::Int32; ny::Int32; nz::Int32
    starts::CuDeviceVector{Int32}
    counts::CuDeviceVector{Int32}
    pts::CuDeviceVector{Int32}        # indices into X/Y/Z (1-based)
    cell_of_point::CuDeviceVector{Int32}
end

struct GPUGridPeriodic
    Lx::Float64; Ly::Float64; Lz::Float64
    hx::Float64; hy::Float64; hz::Float64
    nx::Int32; ny::Int32; nz::Int32
    starts::CuDeviceVector{Int32}
    counts::CuDeviceVector{Int32}
    pts::CuDeviceVector{Int32}
    cell_of_point::CuDeviceVector{Int32}
end

# Build GPU mirrors from your existing CPU grids
function to_gpu_grid(G)::GPUGridNonPeriodic
    # Build a point->cell map once (host) and upload
    cell_of_point = Vector{Int32}(undef, length(G.pts))
    # inverse mapping: we can fill by scanning all cells; but we know G.pts lists points once
    # To find each point's cell, we loop cells and assign
    for cid in 1:length(G.counts)
        s = G.starts[cid]; c = G.counts[cid]
        @inbounds for t in s:(s+c-1)
            p = G.pts[t]
            cell_of_point[p] = Int32(cid)
        end
    end
    GPUGridNonPeriodic(
        G.x0, G.y0, G.z0, G.invcs,
        Int32(G.nx), Int32(G.ny), Int32(G.nz),
        CuArray(Int32.(G.starts)),
        CuArray(Int32.(G.counts)),
        CuArray(Int32.(G.pts)),
        CuArray(cell_of_point),
    )
end

function to_gpu_grid(G; Lx::Float64, Ly::Float64, Lz::Float64)::GPUGridPeriodic
    cell_of_point = Vector{Int32}(undef, length(G.pts))
    for cid in 1:length(G.counts)
        s = G.starts[cid]; c = G.counts[cid]
        @inbounds for t in s:(s+c-1)
            p = G.pts[t]
            cell_of_point[p] = Int32(cid)
        end
    end
    GPUGridPeriodic(
        Lx, Ly, Lz, 0.5Lx, 0.5Ly, 0.5Lz,
        Int32(G.nx), Int32(G.ny), Int32(G.nz),
        CuArray(Int32.(G.starts)),
        CuArray(Int32.(G.counts)),
        CuArray(Int32.(G.pts)),
        CuArray(cell_of_point),
    )
end

# Decode (ix,iy,iz) from linear cell id (1-based)
@inline function gpu_decode_cell(cid::Int32, ny::Int32, nz::Int32)
    t = cid - 1
    ix = t ÷ (ny*nz)
    r  = t % (ny*nz)
    iy = r ÷ nz
    iz = r % nz
    return ix, iy, iz
end

# Clamp neighbor ranges (non-periodic)
@inline function clamp_range(a0::Int32, a1::Int32, n::Int32)
    lo = a0 < 0 ? Int32(0) : a0
    hi = a1 >= n ? n-1 : a1
    return lo, hi
end

# Periodic wrap
@inline function wrap(a::Int32, n::Int32)
    v = a % n
    return v < 0 ? v + n : v
end

# Min-image for one component
@inline function minimg_d(Δ::Float64, halfL::Float64, L::Float64)
    if Δ >  halfL
        Δ - L
    elseif Δ <= -halfL
        Δ + L
    else
        Δ
    end
end

# ------------------ Kernels ------------------

const MAX_NEI = 2048  # per-point neighbor buffer cap

# Non-periodic, true LOS, GPU
function _kernel_nonperiodic!(
    X::CuDeviceVector{Float64}, Y::CuDeviceVector{Float64}, Z::CuDeviceVector{Float64},
    G::GPUGridNonPeriodic,
    rmin2::Float64, rmax2::Float64, μmax2::Float64,
    rmin::Float64, invΔr::Float64, Nr::Int32,
    invΔμ::Float64, Nμ::Int32,
    H::CuDeviceArray{Int32,4},
)
    i = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    Ni = Int32(length(X))
    if i < 1 || i > Ni
        return
    end

    xi = X[i]; yi = Y[i]; zi = Z[i]

    # Decode this point's cell & neighbor range
    cid = G.cell_of_point[i]
    ix, iy, iz = gpu_decode_cell(cid, G.ny, G.nz)
    x0, x1 = clamp_range(ix-1, ix+1, G.nx)
    y0, y1 = clamp_range(iy-1, iy+1, G.ny)
    z0, z1 = clamp_range(iz-1, iz+1, G.nz)

    # local neighbor buffer (indices j)
    nei = MVector{MAX_NEI, Int32}(0)
    m = Int32(0)

    # gather neighbors j>i with (r, μ) cuts
    @inbounds for cx in x0:x1, cy in y0:y1, cz in z0:z1
        ccid = (cx * G.ny) * G.nz + cy * G.nz + cz + 1
        s = G.starts[ccid]; c = G.counts[ccid]
        for t in s:(s+c-1)
            j = G.pts[t]
            if j > i
                dx = Float64(X[j] - xi)
                dy = Float64(Y[j] - yi)
                dz = Float64(Z[j] - zi)
                r2 = dx*dx + dy*dy + dz*dz
                if r2 >= rmin2 && r2 < rmax2
                    num2, den2 = gpu_mu_true_los_numden2(xi,yi,zi, Float64(X[j]),Float64(Y[j]),Float64(Z[j]))
                    if (den2 == 0.0) || (num2 < μmax2 * den2)
                        if m < MAX_NEI
                            m += 1
                            Base.setindex!(nei, j, m) # write into NTuple
                        end
                    end
                end
            end
        end
    end

    # Triangles from neighbor list
    if m >= 2
        for a in Int32(1):(m-1)
            j = nei[a]
            xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
            dxij = xj - xi; dyij = yj - yi; dzij = zj - zi
            r2_ij = dxij*dxij + dyij*dyij + dzij*dzij
            mdij_n, mdij_d = gpu_mu_true_los_numden2(xi,yi,zi, xj,yj,zj)

            for b in (a+1):m
                k = nei[b]
                xk = Float64(X[k]); yk = Float64(Y[k]); zk = Float64(Z[k])

                dxjk = xk - xj; dyjk = yk - yj; dzjk = zk - zj
                r2_jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk
                if r2_jk < rmin2 || r2_jk >= rmax2
                    continue
                end
                mdjk_n, mdjk_d = gpu_mu_true_los_numden2(xj,yj,zj, xk,yk,zk)
                if !(mdjk_d == 0.0 || mdjk_n < μmax2 * mdjk_d)
                    continue
                end

                dxik = xi - xk; dyik = yi - yk; dzik = zi - zk
                r2_ik = dxik*dxik + dyik*dyik + dzik*dzik
                mdik_n, mdik_d = gpu_mu_true_los_numden2(xi,yi,zi, xk,yk,zk)

                r_a, r_b, a_n, a_d, b_n, b_d =
                    gpu_two_shortest_using_r2(r2_ij, r2_jk, r2_ik,
                                              mdij_n, mdij_d,
                                              mdjk_n, mdjk_d,
                                              mdik_n, mdik_d)
                μ_a = gpu_mu_from_numden2(a_n, a_d)
                μ_b = gpu_mu_from_numden2(b_n, b_d)

                ir = gpu_bin_index_inv(r_a, rmin, invΔr, Nr); ir == 0 && continue
                jr = gpu_bin_index_inv(r_b, rmin, invΔr, Nr); jr == 0 && continue
                im = gpu_bin_index_inv(μ_a, 0.0,  invΔμ, Nμ); im == 0 && continue
                jm = gpu_bin_index_inv(μ_b, 0.0,  invΔμ, Nμ); jm == 0 && continue
                @inbounds CUDA.@atomic H[ir, jr, im, jm] += Int32(1)
            end
        end
    end
    return
end

# Periodic, z-LOS, GPU
function _kernel_periodic!(
    X::CuDeviceVector{Float64}, Y::CuDeviceVector{Float64}, Z::CuDeviceVector{Float64},
    G::GPUGridPeriodic,
    rmin2::Float64, rmax2::Float64, μmax2::Float64,
    rmin::Float64, invΔr::Float64, Nr::Int32,
    invΔμ::Float64, Nμ::Int32,
    H::CuDeviceArray{Int32,4},
)
    i = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    Ni = Int32(length(X))
    if i < 1 || i > Ni
        return
    end

    xi = X[i]; yi = Y[i]; zi = Z[i]

    cid = G.cell_of_point[i]
    ix, iy, iz = gpu_decode_cell(cid, G.ny, G.nz)

    # neighbor shells with periodic wrap
    # local neighbor buffer
    nei = MVector{MAX_NEI, Int32}(0)
    m = Int32(0)

    for dx in Int32(-1):Int32(1), dy in Int32(-1):Int32(1), dz in Int32(-1):Int32(1)
        jx = wrap(ix + dx, G.nx)
        jy = wrap(iy + dy, G.ny)
        jz = wrap(iz + dz, G.nz)
        ccid = (jx * G.ny) * G.nz + jy * G.nz + jz + 1
        s = G.starts[ccid]; c = G.counts[ccid]
        @inbounds for t in s:(s+c-1)
            j = G.pts[t]
            if j > i
                dxij = Float64(X[j] - xi)
                dyij = Float64(Y[j] - yi)
                dzij = Float64(Z[j] - zi)
                # min-image
                dxij = minimg_d(dxij, G.hx, G.Lx)
                dyij = minimg_d(dyij, G.hy, G.Ly)
                dzij = minimg_d(dzij, G.hz, G.Lz)
                r2 = dxij*dxij + dyij*dyij + dzij*dzij
                if r2 >= rmin2 && r2 < rmax2
                    num2, den2 = gpu_mu_z_los_numden2(dxij, dyij, dzij)
                    if (den2 == 0.0) || (num2 < μmax2 * den2)
                        if m < MAX_NEI
                            m += 1
                            Base.setindex!(nei, j, m)
                        end
                    end
                end
            end
        end
    end

    if m >= 2
        for a in Int32(1):(m-1)
            j = nei[a]
            xj = Float64(X[j]); yj = Float64(Y[j]); zj = Float64(Z[j])
            dxij = minimg_d(xj - xi, G.hx, G.Lx)
            dyij = minimg_d(yj - yi, G.hy, G.Ly)
            dzij = minimg_d(zj - zi, G.hz, G.Lz)
            r2_ij = dxij*dxij + dyij*dyij + dzij*dzij
            mdij_n, mdij_d = gpu_mu_z_los_numden2(dxij, dyij, dzij)

            for b in (a+1):m
                k = nei[b]
                xk = Float64(X[k]); yk = Float64(Y[k]); zk = Float64(Z[k])

                dxjk = minimg_d(xk - xj, G.hx, G.Lx)
                dyjk = minimg_d(yk - yj, G.hy, G.Ly)
                dzjk = minimg_d(zk - zj, G.hz, G.Lz)
                r2_jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk
                if r2_jk < rmin2 || r2_jk >= rmax2
                    continue
                end
                mdjk_n, mdjk_d = gpu_mu_z_los_numden2(dxjk, dyjk, dzjk)
                if !(mdjk_d == 0.0 || mdjk_n < μmax2 * mdjk_d)
                    continue
                end

                dxik = minimg_d(xi - xk, G.hx, G.Lx)
                dyik = minimg_d(yi - yk, G.hy, G.Ly)
                dzik = minimg_d(zi - zk, G.hz, G.Lz)
                r2_ik = dxik*dxik + dyik*dyik + dzik*dzik
                mdik_n, mdik_d = gpu_mu_z_los_numden2(dxik, dyik, dzik)

                r_a, r_b, a_n, a_d, b_n, b_d =
                    gpu_two_shortest_using_r2(r2_ij, r2_jk, r2_ik,
                                              mdij_n, mdij_d,
                                              mdjk_n, mdjk_d,
                                              mdik_n, mdik_d)
                μ_a = gpu_mu_from_numden2(a_n, a_d)
                μ_b = gpu_mu_from_numden2(b_n, b_d)

                ir = gpu_bin_index_inv(r_a, rmin, invΔr, Nr); ir == 0 && continue
                jr = gpu_bin_index_inv(r_b, rmin, invΔr, Nr); jr == 0 && continue
                im = gpu_bin_index_inv(μ_a, 0.0,  invΔμ, Nμ); im == 0 && continue
                jm = gpu_bin_index_inv(μ_b, 0.0,  invΔμ, Nμ); jm == 0 && continue
                @inbounds CUDA.@atomic H[ir, jr, im, jm] += Int32(1)
            end
        end
    end
    return
end

# ------------------ Public GPU entry points ------------------

"""
    count_triangles_grid_gpu!(X,Y,Z; rmin, rmax, Nr, μmax, Nμ=2, cellsize=rmax, threads=256)

Non-periodic, true LOS. Builds the CPU grid, mirrors it to GPU, and runs a CUDA kernel.
Returns `HistR12R23Mu12Mu13{Int32}` on the host.
"""
function count_triangles_grid_gpu!(X::AbstractVector{<:Real},
                                   Y::AbstractVector{<:Real},
                                   Z::AbstractVector{<:Real};
                                   rmin::Real, rmax::Real, Nr::Integer,
                                   μmax::Real, Nμ::Integer=2,
                                   cellsize::Real=rmax,
                                   threads::Integer=256)

    N = length(X); @assert length(Y)==N && length(Z)==N
    # constants
    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf*μmaxf
    NrI = Int32(Nr); NμI = Int32(Nμ)
    invΔr = Float64(Nr) / (rmaxf - rminf)
    invΔμ = Float64(Nμ) / μmaxf
    r2min = rminf*rminf
    r2max = rmaxf*rmaxf

    # Build grid on CPU (reusing your function) and mirror
    Gcpu = TriCo.build_grid_nonperiodic(X, Y, Z, cellsize)
    Ggpu = to_gpu_grid(Gcpu)

    # Upload coordinates
    dX = CuArray(Float64.(X)); dY = CuArray(Float64.(Y)); dZ = CuArray(Float64.(Z))

    # Device histogram
    dH = CuArray(zeros(Int32, Int(Nr), Int(Nr), Int(Nμ), Int(Nμ)))

    blocks = cld(N, threads)
    @cuda threads=threads blocks=blocks _kernel_nonperiodic!(dX, dY, dZ, Ggpu,
        r2min, r2max, μmax2,
        rminf, invΔr, NrI,
        invΔμ, NμI,
        dH)

    Hh = Array(dH)
    return HistR12R23Mu12Mu13(Int(Nr), Int(Nμ); T=Int32) |> (H -> (H.h .= Hh; H))
end

"""
    count_triangles_periodic_grid_gpu!(X,Y,Z; Lx,Ly,Lz, rmin,rmax,Nr, μmax, Nμ=2, cellsize=rmax, threads=256)

Periodic, min-image (z is LOS). CPU grid build, GPU enumeration and binning.
"""
function count_triangles_periodic_grid_gpu!(X::AbstractVector{<:Real},
                                            Y::AbstractVector{<:Real},
                                            Z::AbstractVector{<:Real};
                                            Lx::Real, Ly::Real, Lz::Real,
                                            rmin::Real, rmax::Real, Nr::Integer,
                                            μmax::Real, Nμ::Integer=2,
                                            cellsize::Real=rmax,
                                            threads::Integer=256)

    N = length(X); @assert length(Y)==N && length(Z)==N
    rminf = Float64(rmin); rmaxf = Float64(rmax)
    μmaxf = Float64(μmax); μmax2 = μmaxf*μmaxf
    NrI = Int32(Nr); NμI = Int32(Nμ)
    invΔr = Float64(Nr) / (rmaxf - rminf)
    invΔμ = Float64(Nμ) / μmaxf
    r2min = rminf*rminf
    r2max = rmaxf*rmaxf

    # Build periodic grid on CPU then mirror
    Gcpu = TriCo.build_grid_periodic(X, Y, Z; Lx=Lx, Ly=Ly, Lz=Lz, cellsize=cellsize)
    Ggpu = to_gpu_grid(Gcpu; Lx=Float64(Lx), Ly=Float64(Ly), Lz=Float64(Lz))

    dX = CuArray(Float64.(X)); dY = CuArray(Float64.(Y)); dZ = CuArray(Float64.(Z))
    dH = CuArray(zeros(Int32, Int(Nr), Int(Nr), Int(Nμ), Int(Nμ)))

    blocks = cld(N, threads)
    @cuda threads=threads blocks=blocks _kernel_periodic!(dX, dY, dZ, Ggpu,
        r2min, r2max, μmax2,
        rminf, invΔr, NrI,
        invΔμ, NμI,
        dH)

    Hh = Array(dH)
    return HistR12R23Mu12Mu13(Int(Nr), Int(Nμ); T=Int32) |> (H -> (H.h .= Hh; H))
end

end # module

