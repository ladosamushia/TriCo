#!/usr/bin/env julia
# --------------------------------------------------------------
# Use TriCo.read_xyz_fits + count_triangles!(...) or
# count_triangles_periodic!(...) to build histogram from FITS.
# --------------------------------------------------------------

using TriCo
using Base.Threads

# --- tiny arg parser (same as before) ---
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

function main()
    d = parse_args(ARGS)
    if !haskey(d, "--fits")
        println("Usage: julia --project=. scripts/fits_to_hist.jl --fits FILE.fits [options]")
        return
    end

    fits   = d["--fits"]
    xcol   = getopt(d, "--xcol", "X")
    ycol   = getopt(d, "--ycol", "Y")
    zcol   = getopt(d, "--zcol", "Z")
    hdu    = getopt_int(d, "--hdu", 2)

    rmin   = getopt_float(d, "--rmin", 5.0)
    rmax   = getopt_float(d, "--rmax", 60.0)
    Nr     = getopt_int(d, "--Nr", 55)
    mumax  = getopt_float(d, "--mumax", 0.9)
    Nmu    = getopt_int(d, "--Nmu", 2)
    cellsize = getopt_float(d, "--cellsize", rmax)

    periodic = getopt_bool(d, "--periodic", false)
    Lx = getopt_float(d, "--Lx", 2000.0)
    Ly = getopt_float(d, "--Ly", 2000.0)
    Lz = getopt_float(d, "--Lz", 2000.0)

    println("Threads: ", Threads.nthreads())
    println("Reading FITS: $fits (HDU=$hdu, cols: $xcol,$ycol,$zcol)")
    X, Y, Z = read_xyz_fits(fits; xcol=xcol, ycol=ycol, zcol=zcol, hdu=hdu)
    println("Loaded ", length(X), " rows.")

    if periodic
        X .= mod.(X, Lx); Y .= mod.(Y, Ly); Z .= mod.(Z, Lz)
        H = count_triangles_periodic!(X, Y, Z;
                                      Lx=Lx, Ly=Ly, Lz=Lz,
                                      rmin=rmin, rmax=rmax, Nr=Nr,
                                      μmax=mumax, Nμ=Nmu, cellsize=cellsize)
    else
        H = count_triangles!(X, Y, Z;
                             rmin=rmin, rmax=rmax, Nr=Nr,
                             μmax=mumax, Nμ=Nmu, cellsize=cellsize)
    end

    println("Histogram size : ", size(H.h))
    println("Total triangles: ", sum(H.h))
    println("Nonzero bins   : ", count(!=(0), H.h), " / ", Nr*Nr*Nmu*Nmu)
end

main()

