#!/usr/bin/env julia
using Random
using TriCo
using Base.Threads

# ---------- Parameters ----------
N    = parse(Int,    get(ENV, "TRICO_N",    "100000"))
L    = parse(Float64,get(ENV, "TRICO_L",    "2000"))
rmin = parse(Float64,get(ENV, "TRICO_RMIN", "5.0"))
rmax = parse(Float64,get(ENV, "TRICO_RMAX", "60.0"))
Nr   = parse(Int,    get(ENV, "TRICO_NR",   "55"))
μmax = parse(Float64,get(ENV, "TRICO_MUMAX","0.9"))   # default 0.9
Nμ   = parse(Int,    get(ENV, "TRICO_NMU",  "2"))
CELL = parse(Float64,get(ENV, "TRICO_CELL", string(rmax)))
SEED = parse(Int,    get(ENV, "TRICO_SEED", "12345"))

Random.seed!(SEED)

println("Threads: ", Threads.nthreads())
println("Generating $N random points in cube of size $L...")

# Draw on (0,1) and scale by L
X = L .* rand(N)
Y = L .* rand(N)
Z = L .* rand(N)

println("Running count_triangles! ...")
t = @elapsed begin
    H4 = count_triangles!(X, Y, Z; rmin=rmin, rmax=rmax, Nr=Nr, μmax=μmax, Nμ=Nμ, cellsize=CELL)
    total = sum(H4.h)
    nz = count(!=(0), H4.h)
    println("Histogram size   : ", size(H4.h))          # (Nr, Nr, Nμ, Nμ)
    println("Nonzero bins     : ", nz, " / ", Nr*Nr*Nμ*Nμ)
    println("Total triangles  : ", total)
end
println("Elapsed time     : $(round(t, digits=3)) s")
println("Done.")

