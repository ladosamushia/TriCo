#!/usr/bin/env julia
# scripts/fits_to_hist.jl  (single + mixed catalog support)
# --------------------------------------------------------------
# Build a histogram from FITS table(s) using grid-accelerated kernels:
#   - TriCo.count_triangles_grid!                  (non-periodic, true LOS)
#   - TriCo.count_triangles_periodic_grid!         (periodic, min-image, z-LOS)
#   - TriCo.count_triangles_mixed!(A,B,C)          (non-periodic, mixed catalogs)
#   - TriCo.count_triangles_mixed_periodic!(A,B,C) (periodic, mixed catalogs)
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
    Usage (single catalog):
      julia --project=. scripts/fits_to_hist.jl --fits FILE.fits [options]

    Usage (mixed catalogs):
      julia --project=. scripts/fits_to_hist.jl --fitsA A.fits [--fitsB B.fits] [--fitsC C.fits] [--pattern AAB|ABB|ABC] [options]

      If two catalogs are given (A and B) and --pattern is not provided, defaults to AAB (two from A, one from B).
      If three catalogs are given (A, B, C) and --pattern is not provided, defaults to ABC.

    Required:
      --fits PATH                  (single-catalog mode) path to FITS file
      --fitsA PATH                 (mixed mode) path to FITS for catalog A
      --fitsB PATH                 (optional) catalog B
      --fitsC PATH                 (optional) catalog C

    Optional (FITS columns; apply to all unless A/B/C suffix provided):
      --xcol X                     column name for X (default: X)        [or --xcolA, --xcolB, --xcolC]
      --ycol Y                     column name for Y (default: Y)        [or --ycolA, --ycolB, --ycolC]
      --zcol Z                     column name for Z (default: Z)        [or --zcolA, --zcolB, --zcolC]
      --hdu  2                     HDU number (default: 2)               [or --hduA, --hduB, --hduC]

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

    Mixed composition:
      --pattern AAB|ABB|ABC        which vertices come from which catalogs (default depends on inputs)

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
        (hdu < 1 || hdu > length(f)) && error("HDU=$(hdu) out of range for $(path)")
        t = f[hdu]
        X = read(t, xcol)
        Y = read(t, ycol)
        Z = read(t, zcol)
    end
    return collect(Float64.(X)), collect(Float64.(Y)), collect(Float64.(Z))
end

# Helper to get possibly A/B/C-suffixed option
function getopt_abc(d::Dict{String,String}, base::String, suf::String, default)
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
    mixed  = any(haskey(d, k) for k in ("--fitsA","--fitsB","--fitsC"))
    if !(single || mixed)
        usage(); return
    end
    if single && mixed
        error("Specify either --fits (single) OR --fitsA/B/C (mixed), not both.")
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

        # compute histogram (GRID kernels)
        println("Computing histogram…")
        H = if periodic
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

        if !isempty(outpath)
            NPZ.npzwrite(outpath, Dict("hist" => H.h))
            println("Saved NPZ (hist only) -> $(outpath)")
        end
        println("Done.")
        return
    end

    # Mixed mode
    fitsA = getopt(d, "--fitsA", "")
    fitsB = getopt(d, "--fitsB", "")
    fitsC = getopt(d, "--fitsC", "")

    nprov = count(!isempty, (fitsA, fitsB, fitsC))
    if nprov < 2
        error("Mixed mode requires at least two catalogs: provide --fitsA and --fitsB (and optionally --fitsC)." )
    end

    pattern = uppercase(getopt(d, "--pattern", nprov == 2 ? "AAB" : "ABC"))
    if !(pattern in ("AAB","ABB","ABC"))
        error("--pattern must be one of AAB, ABB, ABC")
    end
    if pattern == "ABC" && nprov < 3
        error("--pattern=ABC requires --fitsA, --fitsB, and --fitsC.")
    end

    # Columns/HDUs (shared defaults with optional suffix override)
    xcol = getopt(d, "--xcol", "X");    xcolA = getopt_abc(d, "--xcol", "A", xcol); xcolB = getopt_abc(d, "--xcol", "B", xcol); xcolC = getopt_abc(d, "--xcol", "C", xcol)
    ycol = getopt(d, "--ycol", "Y");    ycolA = getopt_abc(d, "--ycol", "A", ycol); ycolB = getopt_abc(d, "--ycol", "B", ycol); ycolC = getopt_abc(d, "--ycol", "C", ycol)
    zcol = getopt(d, "--zcol", "Z");    zcolA = getopt_abc(d, "--zcol", "A", zcol); zcolB = getopt_abc(d, "--zcol", "B", zcol); zcolC = getopt_abc(d, "--zcol", "C", zcol)
    hdu  = getopt_int(d, "--hdu", 2);    hduA  = haskey(d,"--hduA") ? getopt_int(d,"--hduA",2) : hdu
                                          hduB  = haskey(d,"--hduB") ? getopt_int(d,"--hduB",2) : hdu
                                          hduC  = haskey(d,"--hduC") ? getopt_int(d,"--hduC",2) : hdu

    # Read catalogs that are provided
    XA=YA=ZA=XB=YB=ZB=XC=YC=ZC = nothing

    if !isempty(fitsA)
        println("Reading FITS A: $(fitsA) (HDU=$(hduA), cols: $(xcolA),$(ycolA),$(zcolA))")
        XA, YA, ZA = read_xyz_fits(fitsA; xcol=xcolA, ycol=ycolA, zcol=zcolA, hdu=hduA)
        println("Loaded A: ", length(XA))
    end
    if !isempty(fitsB)
        println("Reading FITS B: $(fitsB) (HDU=$(hduB), cols: $(xcolB),$(ycolB),$(zcolB))")
        XB, YB, ZB = read_xyz_fits(fitsB; xcol=xcolB, ycol=ycolB, zcol=zcolB, hdu=hduB)
        println("Loaded B: ", length(XB))
    end
    if !isempty(fitsC)
        println("Reading FITS C: $(fitsC) (HDU=$(hduC), cols: $(xcolC),$(ycolC),$(zcolC))")
        XC, YC, ZC = read_xyz_fits(fitsC; xcol=xcolC, ycol=ycolC, zcol=zcolC, hdu=hduC)
        println("Loaded C: ", length(XC))
    end

    # Build TriCat wrappers according to provided catalogs
    A = TriCo.TriCat(XA,YA,ZA)
    B = TriCo.TriCat(isnothing(XB) ? XA : XB, isnothing(YB) ? YA : YB, isnothing(ZB) ? ZA : ZB)  # fallback to A if B missing
    C = TriCo.TriCat(isnothing(XC) ? (pattern=="ABC" ? error("C required for ABC") : XA) : XC,
                     isnothing(YC) ? (pattern=="ABC" ? error("C required for ABC") : YA) : YC,
                     isnothing(ZC) ? (pattern=="ABC" ? error("C required for ABC") : ZA) : ZC)

    # Periodic wrap if needed
    if periodic
        @inbounds begin
            XA .= mod.(XA, Lx); YA .= mod.(YA, Ly); ZA .= mod.(ZA, Lz)
            if !isnothing(XB); XB .= mod.(XB, Lx); YB .= mod.(YB, Ly); ZB .= mod.(ZB, Lz); end
            if !isnothing(XC); XC .= mod.(XC, Lx); YC .= mod.(YC, Ly); ZC .= mod.(ZC, Lz); end
        end
    end

    println("Computing histogram (mixed; pattern=$(pattern))…")
    H = if periodic
        if pattern == "AAB"
            TriCo.count_triangles_mixed_periodic!(A,A,B; Lx=Lx, Ly=Ly, Lz=Lz,
                rmin=rmin, rmax=rmax, Nr=Nr, μmax=mumax, Nμ=Nmu, cellsize=cellsize)
        elseif pattern == "ABB"
            TriCo.count_triangles_mixed_periodic!(A,B,B; Lx=Lx, Ly=Ly, Lz=Lz,
                rmin=rmin, rmax=rmax, Nr=Nr, μmax=mumax, Nμ=Nmu, cellsize=cellsize)
        else # ABC
            TriCo.count_triangles_mixed_periodic!(A,B,C; Lx=Lx, Ly=Ly, Lz=Lz,
                rmin=rmin, rmax=rmax, Nr=Nr, μmax=mumax, Nμ=Nmu, cellsize=cellsize)
        end
    else
        if pattern == "AAB"
            TriCo.count_triangles_mixed!(A,A,B; rmin=rmin, rmax=rmax, Nr=Nr, μmax=mumax, Nμ=Nmu, cellsize=cellsize)
        elseif pattern == "ABB"
            TriCo.count_triangles_mixed!(A,B,B; rmin=rmin, rmax=rmax, Nr=Nr, μmax=mumax, Nμ=Nmu, cellsize=cellsize)
        else
            TriCo.count_triangles_mixed!(A,B,C; rmin=rmin, rmax=rmax, Nr=Nr, μmax=mumax, Nμ=Nmu, cellsize=cellsize)
        end
    end

    println("Histogram size : ", size(H.h))
    println("Total triangles: ", sum(H.h))
    println("Nonzero bins   : ", count(!=(0), H.h), " / ", Nr*Nr*Nmu*Nmu)

    if !isempty(outpath)
        NPZ.npzwrite(outpath, Dict("hist" => H.h))
        println("Saved NPZ (hist only) -> $(outpath)")
    end

    println("Done.")
end

main()

