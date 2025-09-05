#!/usr/bin/env julia
# scripts/fits_to_hist.jl
# --------------------------------------------------------------
# Build a histogram from a FITS table using grid-accelerated kernels:
#   - TriCo.count_triangles_grid!          (non-periodic, true LOS)
#   - TriCo.count_triangles_periodic_grid! (periodic, min-image, z-LOS)
# Save ONLY the 4D histogram array to .npz (key = "hist").
# --------------------------------------------------------------

using TriCo
using FITSIO
using NPZ
using Base.Threads

# --- tiny arg parser (no deps) ---
function parse_args(args::Vector{String})
    d = Dict{String,String}()
    i = 1
    while i <= length(args)
        a = args[i]
        if startswith(a, "--")
            key = a
            if i < length(args) && !startswith(args[i+1], "--")
                d[key] = args[i+1]; i += 2
            else
                d[key] = "true"; i += 1
            end
        else
            i += 1
        end
    end
    return d
end
getopt(d,k,default) = haskey(d,k) ? d[k] : default
getopt_bool(d,k,default::Bool=false) = lowercase(getopt(d,k,string(default))) in ("1","true","yes","on")
getopt_float(d,k,default) = parse(Float64, getopt(d,k,string(default)))
getopt_int(d,k,default) = parse(Int, getopt(d,k,string(default)))

function usage()
    println("""
    Usage:
      julia --project=. scripts/fits_to_hist.jl --fits FILE.fits [options]

    Required:
      --fits PATH                  path to FITS file

    Optional (FITS):
      --xcol X                     column name for X (default: X)
      --ycol Y                     column name for Y (default: Y)
      --zcol Z                     column name for Z (default: Z)
      --hdu  2                     HDU number of the binary table (default: 2)

    Selection & binning:
      --rmin 5.0                   min r (default 5.0)
      --rmax 60.0                  max r (default 60.0)
      --Nr   55                    number of r bins (default 55)
      --mumax 0.9                  μ_max (default 0.9)
      --Nmu  2                     μ bins (default 2)

    Neighbor grid:
      --cellsize 60.0              grid cell size (default rmax)

    Periodic box (z is LOS):
      --periodic                   enable periodic mode (default OFF)
      --Lx 2000 --Ly 2000 --Lz 2000  box sizes (required if --periodic)

    Output:
      --out triangles.npz          save ONLY the histogram to NPZ (key: "hist")
    """)
end

# --- FITS reader (table HDU, named columns) ---
function read_xyz_fits(path::AbstractString; xcol::AbstractString="X",
                       ycol::AbstractString="Y", zcol::AbstractString="Z",
                       hdu::Integer=2)
    X = Y = Z = nothing
    FITS(path, "r") do f
        (hdu < 1 || hdu > length(f)) && error("HDU=$hdu out of range for $path")
        t = f[hdu]
        X = read(t, xcol)
        Y = read(t, ycol)
        Z = read(t, zcol)
    end
    return collect(Float64.(X)), collect(Float64.(Y)), collect(Float64.(Z))
end

function main()
    d = parse_args(ARGS)
    if !haskey(d, "--fits")
        usage(); return
    end

    # FITS params
    fits   = d["--fits"]
    xcol   = getopt(d, "--xcol", "X")
    ycol   = getopt(d, "--ycol", "Y")
    zcol   = getopt(d, "--zcol", "Z")
    hdu    = getopt_int(d, "--hdu", 2)

    # selection & binning
    rmin   = getopt_float(d, "--rmin", 5.0)
    rmax   = getopt_float(d, "--rmax", 60.0)
    Nr     = getopt_int(d, "--Nr", 55)
    mumax  = getopt_float(d, "--mumax", 0.9)
    Nmu    = getopt_int(d, "--Nmu", 2)

    # grid
    cellsize = getopt_float(d, "--cellsize", rmax)

    # periodic?
    periodic = getopt_bool(d, "--periodic", false)
    Lx = getopt_float(d, "--Lx", 2000.0)
    Ly = getopt_float(d, "--Ly", 2000.0)
    Lz = getopt_float(d, "--Lz", 2000.0)

    # output path (optional)
    outpath = getopt(d, "--out", "")

    # basic sanity
    @assert rmin < rmax "Require rmin < rmax"
    @assert 0.0 < mumax <= 1.0 "μmax must be in (0,1]"
    @assert cellsize >= rmax "cellsize must be ≥ rmax to avoid missing neighbors"
    if periodic
        @assert all(!isnan.((Lx, Ly, Lz))) "Periodic mode requires --Lx, --Ly, --Lz"
    end
    if !isempty(outpath)
        if !endswith(lowercase(outpath), ".npz")
            error("--out must end with .npz")
        end
    end

    println("Threads: ", Threads.nthreads())
    println("Reading FITS: $fits (HDU=$hdu, cols: $xcol,$ycol,$zcol)")
    X, Y, Z = read_xyz_fits(fits; xcol=xcol, ycol=ycol, zcol=zcol, hdu=hdu)
    println("Loaded ", length(X), " rows.")

    # compute histogram (GRID kernels)
    println("Computing histogram…")
    H = if periodic
        # ensure coordinates are within [0,L)
        @inbounds begin
            X .= mod.(X, Lx); Y .= mod.(Y, Ly); Z .= mod.(Z, Lz)
        end
        TriCo.count_triangles_periodic_grid!(X, Y, Z;
            Lx=Lx, Ly=Ly, Lz=Lz,
            rmin=rmin, rmax=rmax, Nr=Nr,
            μmax=mumax, Nμ=Nmu, cellsize=cellsize)
    else
        TriCo.count_triangles_grid!(X, Y, Z;
            rmin=rmin, rmax=rmax, Nr=Nr,
            μmax=mumax, Nμ=Nmu, cellsize=cellsize)
    end

    # summary
    println("Histogram size : ", size(H.h))                # (Nr, Nr, Nμ, Nμ)
    println("Total triangles: ", sum(H.h))
    println("Nonzero bins   : ", count(!=(0), H.h), " / ", Nr*Nr*Nmu*Nmu)

    # save ONLY the histogram array when requested
    if !isempty(outpath)
        NPZ.npzwrite(outpath, Dict("hist" => H.h))
        println("Saved NPZ (hist only) -> $outpath")
    end

    println("Done.")
end

main()

