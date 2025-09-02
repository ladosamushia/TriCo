#!/usr/bin/env julia
# --------------------------------------------------------------
# Periodic random-cube driver for TriCo.count_triangles_periodic_zlos!
# - Box can have different lengths Lx, Ly, Lz
# - z is the line-of-sight: μ^2 = Δz^2 / r^2 (minimum image)
# - All pairs must satisfy rmin < r_ij < rmax and μ_ij < μmax
# - Produces a 4D joint histogram over (r12, r23, μ12, μ13)
# --------------------------------------------------------------

using Random
using Base.Threads
using TriCo   # make sure your module exports count_triangles_periodic!

# ---------------- Parameters (override via ENV) ----------------
N    = parse(Int,    get(ENV, "TRICO_N",    "200000"))

Lx   = parse(Float64,get(ENV, "TRICO_LX",   "2000"))
Ly   = parse(Float64,get(ENV, "TRICO_LY",   "2000"))
Lz   = parse(Float64,get(ENV, "TRICO_LZ",   "2000"))

rmin = parse(Float64,get(ENV, "TRICO_RMIN", "5.0"))
rmax = parse(Float64,get(ENV, "TRICO_RMAX", "60.0"))
Nr   = parse(Int,    get(ENV, "TRICO_NR",   "55"))

μmax = parse(Float64,get(ENV, "TRICO_MUMAX","0.9"))     # default 0.9, as requested
Nμ   = parse(Int,    get(ENV, "TRICO_NMU",  "2"))       # default 2 μ bins

CELL = parse(Float64,get(ENV, "TRICO_CELL", string(rmax)))  # grid cell size ~ rmax

SEED = parse(Int,    get(ENV, "TRICO_SEED", "12345"))
Random.seed!(SEED)

println("Threads: ", Threads.nthreads())
println("Generating $N random points in periodic box (Lx=$Lx, Ly=$Ly, Lz=$Lz)...")

# Points uniform in [0, L) per axis
X = Lx .* rand(N)
Y = Ly .* rand(N)
Z = Lz .* rand(N)

println("Running count_triangles_periodic! ...")
t = @elapsed begin
    H4 = count_triangles_periodic!(X, Y, Z;
                                        Lx=Lx, Ly=Ly, Lz=Lz,
                                        rmin=rmin, rmax=rmax, Nr=Nr,
                                        μmax=μmax, Nμ=Nμ, cellsize=CELL)

    total = sum(H4.h)
    nz = count(!=(0), H4.h)

    println("Histogram size   : ", size(H4.h))          # (Nr, Nr, Nμ, Nμ)
    println("Nonzero bins     : ", nz, " / ", Nr*Nr*Nμ*Nμ)
    println("Total triangles  : ", total)
end
println("Elapsed time     : $(round(t, digits=3)) s")
println("Done.")

