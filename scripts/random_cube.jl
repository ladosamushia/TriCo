#!/usr/bin/env julia
# scripts/random_cube.jl
#
# Driver for grid-accelerated triangle counting.
# Uses:
#   - count_triangles_grid!                (non-periodic, true LOS)
#   - count_triangles_periodic_grid!       (periodic, min-image, z-LOS)
#
# Environment overrides:
#   TRICO_N, TRICO_L, TRICO_RMIN, TRICO_RMAX, TRICO_NR,
#   TRICO_MUMAX, TRICO_NMU, TRICO_CELL, TRICO_SEED, TRICO_PERIODIC
#
# Examples:
#   julia --project=. scripts/random_cube.jl
#   TRICO_PERIODIC=1 JULIA_NUM_THREADS=8 julia --project=. scripts/random_cube.jl
#   TRICO_N=100000 TRICO_L=1000 TRICO_RMIN=2 TRICO_RMAX=50 TRICO_NR=40 TRICO_MUMAX=0.8 TRICO_NMU=4 julia --project=. scripts/random_cube.jl

using Random
using TriCo
using Base.Threads

# ---------- Parameters ----------
N        = parse(Int,     get(ENV, "TRICO_N",        "200000"))
L        = parse(Float64, get(ENV, "TRICO_L",        "2000"))
rmin     = parse(Float64, get(ENV, "TRICO_RMIN",     "5.0"))
rmax     = parse(Float64, get(ENV, "TRICO_RMAX",     "60.0"))
Nr       = parse(Int,     get(ENV, "TRICO_NR",       "55"))
μmax     = parse(Float64, get(ENV, "TRICO_MUMAX",    "0.9"))   # linear μ threshold
Nμ       = parse(Int,     get(ENV, "TRICO_NMU",      "2"))
CELL     = parse(Float64, get(ENV, "TRICO_CELL",     string(rmax)))
SEED     = parse(Int,     get(ENV, "TRICO_SEED",     "12345"))
PERIODIC = parse(Int,     get(ENV, "TRICO_PERIODIC", "0")) != 0  # 1 => periodic grid

Random.seed!(SEED)

println("Threads          : ", Threads.nthreads())
println("Points           : $N")
println("Box size (L)     : $L")
println("r bins / μ bins  : $Nr / $Nμ")
println("r in [$(rmin), $(rmax)), μ < $(μmax)")
println("Cell size        : $(CELL)  (suggest ~ rmax)")
println("Mode             : ", PERIODIC ? "periodic grid (min-image, z-LOS)" : "non-periodic grid (true LOS)")

# Generate uniform points in [0, L)
X = L .* rand(N)
Y = L .* rand(N)
Z = L .* rand(N)

# ---------- Run ----------
if PERIODIC
    println("Running count_triangles_periodic_grid! ...")
else
    println("Running count_triangles_grid! ...")
end

t = @elapsed begin
    if PERIODIC
        H4 = count_triangles_periodic_grid!(X, Y, Z;
            Lx=L, Ly=L, Lz=L,
            rmin=rmin, rmax=rmax, Nr=Nr,
            μmax=μmax, Nμ=Nμ, cellsize=CELL)
        total = sum(H4.h)
        nz = count(!=(0), H4.h)
        println("Histogram size   : ", size(H4.h))          # (Nr, Nr, Nμ, Nμ)
        println("Nonzero bins     : ", nz, " / ", Nr*Nr*Nμ*Nμ)
        println("Total triangles  : ", total)
    else
        H4 = count_triangles_grid!(X, Y, Z;
            rmin=rmin, rmax=rmax, Nr=Nr,
            μmax=μmax, Nμ=Nμ, cellsize=CELL)
        total = sum(H4.h)
        nz = count(!=(0), H4.h)
        println("Histogram size   : ", size(H4.h))          # (Nr, Nr, Nμ, Nμ)
        println("Nonzero bins     : ", nz, " / ", Nr*Nr*Nμ*Nμ)
        println("Total triangles  : ", total)
    end
end

println("Elapsed time     : $(round(t, digits=3)) s")
tri_per_sec = (t > 0) ? round(total / t, digits=2) : NaN
println("Throughput       : $(tri_per_sec) tris/s (binned)")
println("Done.")

