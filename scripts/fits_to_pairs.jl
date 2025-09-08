#!/usr/bin/env julia
# scripts/fits_to_pairs.jl
# --------------------------------------------------------------
# Build a PAIRS histogram from one or two FITS table(s) using
# grid-accelerated kernels in TriCo:
#   - TriCo.count_pairs_grid!                    (non-periodic, true LOS)
#   - TriCo.count_pairs_periodic_grid!           (periodic, min-image, z-LOS)
#   - TriCo.count_pairs_cross_grid!              (non-periodic, cross catalogs)
#   - TriCo.count_pairs_cross_periodic_grid!     (periodic, cross catalogs)
# Save ONLY the 2D histogram array to .npz (key = "hist").
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
    Usage (single catalog):
      julia --project=. scripts/fits_to_pairs.jl --fits FILE.fits [options]

    Usage (cross catalogs):
      julia --project=. scripts/fits_to_pairs.jl --fitsA A.fits --fitsB B.fits [options]

    Required:
      --fits PATH        (single-catalog mode) path to FITS file
      --fitsA PATH       (cross mode) path to first FITS catalog
      --fitsB PATH       (cross mode) path to second FITS catalog

    Optional (FITS columns; apply to all unless A/B suffix provided):
      --xcol X           column name for X (default: X)         [or --xcolA, --xcolB]
      --ycol Y           column name for Y (default: Y)         [or --ycolA, --ycolB]
      --zcol Z           column name for Z (default: Z)         [or --zcolA, --zcolB]
      --hdu  2           HDU number (default: 2)                [or --hduA, --hduB]

    Selection & binning:
      --rmin 5.0         min r (default 5.0)
      --rmax 60.0        max r (default 60.0)
      --Nr   55          number of r bins (default 55)
      --mumax 0.9        μ_max (default 0.9)
      --Nmu  2           μ bins (default 2)

    Neighbor grid:
      --cellsize 60.0    grid cell size (default rmax)

    Periodic box (z is LOS):
      --periodic         enable periodic mode (default OFF)
      --Lx 2000 --Ly 2000 --Lz 2000  box sizes (required if --periodic)

    Output:
      --out pairs.npz    save ONLY the histogram to NPZ (key: "hist")
    """)
end

# --- FITS reader (table HDU, named columns) ---
function read_xyz_fits(path::AbstractString; xcol::AbstractString="X",
                       ycol::AbstractString="Y", zcol::AbstractString="Z",
                       hdu::Integer=2)
    X = Y = Z = nothing
    FITS(path, "r") do f
        (hdu < 1 || hdu > length(f)) && error("HDU=$(hdu) out of range for $(path)")
        t = f[hdu]
        X = read(t, xcol)
        Y = read(t, ycol)
        Z = read(t, zcol)
    end
    return collect(Float64.(X)), collect(Float64.(Y)), collect(Float64.(Z))
end

# Helper to get possibly A/B-suffixed option
function getopt_ab(d::Dict{String,String}, base::String, suf::String, default)
    key = base * suf
    if haskey(d, key)
        v = d[key]
        if base in ("--hdu",)
            return parse(Int, v)
        else
            return v
        end
    else
        return default
    end
end

function main()
    d = parse_args(ARGS)
    single = haskey(d, "--fits")
    cross  = all(haskey(d, k) for k in ("--fitsA","--fitsB"))
    if !(single || cross)
        usage(); return
    end
    if single && cross
        error("Specify either --fits (single) OR --fitsA/--fitsB (cross), not both.")
    end

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

    if single
        # FITS params (single)
        fits = d["--fits"]
        xcol = getopt(d, "--xcol", "X")
        ycol = getopt(d, "--ycol", "Y")
        zcol = getopt(d, "--zcol", "Z")
        hdu  = getopt_int(d, "--hdu", 2)

        println("Reading FITS: $(fits) (HDU=$(hdu), cols: $(xcol),$(ycol),$(zcol))")
        X, Y, Z = read_xyz_fits(fits; xcol=xcol, ycol=ycol, zcol=zcol, hdu=hdu)
        println("Loaded ", length(X), " rows.")

        # compute pairs histogram
        println("Computing pair histogram…")
        H = if periodic
            @inbounds begin
                X .= mod.(X, Lx); Y .= mod.(Y, Ly); Z .= mod.(Z, Lz)
            end
            TriCo.count_pairs_periodic_grid!(X, Y, Z;
                Lx=Lx, Ly=Ly, Lz=Lz,
                rmin=rmin, rmax=rmax, Nr=Nr,
                μmax=mumax, Nμ=Nmu, cellsize=cellsize)
        else
            TriCo.count_pairs_grid!(X, Y, Z;
                rmin=rmin, rmax=rmax, Nr=Nr,
                μmax=mumax, Nμ=Nmu, cellsize=cellsize)
        end

        # summary
        println("Histogram size : ", size(H.h))      # (Nr, Nμ)
        println("Total pairs    : ", sum(H.h))
        println("Nonzero bins   : ", count(!=(0), H.h), " / ", Nr*Nmu)

        if !isempty(outpath)
            NPZ.npzwrite(outpath, Dict("hist" => H.h))
            println("Saved NPZ (hist only) -> $(outpath)")
        end
        println("Done.")
        return
    end

    # Cross mode
    fitsA = getopt(d, "--fitsA", "")
    fitsB = getopt(d, "--fitsB", "")

    # Columns/HDUs (shared defaults with optional suffix override)
    xcol = getopt(d, "--xcol", "X");    xcolA = getopt_ab(d, "--xcol", "A", xcol); xcolB = getopt_ab(d, "--xcol", "B", xcol)
    ycol = getopt(d, "--ycol", "Y");    ycolA = getopt_ab(d, "--ycol", "A", ycol); ycolB = getopt_ab(d, "--ycol", "B", ycol)
    zcol = getopt(d, "--zcol", "Z");    zcolA = getopt_ab(d, "--zcol", "A", zcol); zcolB = getopt_ab(d, "--zcol", "B", zcol)
    hdu  = getopt_int(d, "--hdu", 2);    hduA  = haskey(d,"--hduA") ? getopt_int(d,"--hduA",2) : hdu
                                         hduB  = haskey(d,"--hduB") ? getopt_int(d,"--hduB",2) : hdu

    println("Reading FITS A: $(fitsA) (HDU=$(hduA), cols: $(xcolA),$(ycolA),$(zcolA))")
    XA, YA, ZA = read_xyz_fits(fitsA; xcol=xcolA, ycol=ycolA, zcol=zcolA, hdu=hduA)
    println("Loaded A: ", length(XA))

    println("Reading FITS B: $(fitsB) (HDU=$(hduB), cols: $(xcolB),$(ycolB),$(zcolB))")
    XB, YB, ZB = read_xyz_fits(fitsB; xcol=xcolB, ycol=ycolB, zcol=zcolB, hdu=hduB)
    println("Loaded B: ", length(XB))

    # Periodic wrap if needed
    if periodic
        @inbounds begin
            XA .= mod.(XA, Lx); YA .= mod.(YA, Ly); ZA .= mod.(ZA, Lz)
            XB .= mod.(XB, Lx); YB .= mod.(YB, Ly); ZB .= mod.(ZB, Lz)
        end
    end

    println("Computing pair histogram (cross)…")
    H = if periodic
        TriCo.count_pairs_cross_periodic_grid!(XA, YA, ZA, XB, YB, ZB;
                                               Lx=Lx, Ly=Ly, Lz=Lz,
                                               rmin=rmin, rmax=rmax, Nr=Nr,
                                               μmax=mumax, Nμ=Nmu, cellsize=cellsize)
    else
        TriCo.count_pairs_cross_grid!(XA, YA, ZA, XB, YB, ZB;
                                      rmin=rmin, rmax=rmax, Nr=Nr,
                                      μmax=mumax, Nμ=Nmu, cellsize=cellsize)
    end

    println("Histogram size : ", size(H.h))   # (Nr, Nμ)
    println("Total pairs    : ", sum(H.h))
    println("Nonzero bins   : ", count(!=(0), H.h), " / ", Nr*Nmu)

    if !isempty(outpath)
        NPZ.npzwrite(outpath, Dict("hist" => H.h))
        println("Saved NPZ (hist only) -> $(outpath)")
    end

    println("Done.")
end

main()

