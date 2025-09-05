# src/io_save.jl

using HDF5

"""
    save_hist_h5(path, H4; rmin, rmax, Nr, mumax, Nmu, periodic=false, Lx=NaN, Ly=NaN, Lz=NaN)

Save the histogram to a single HDF5 file with:
- datasets: "hist" (Nr×Nr×Nμ×Nμ), "r_edges" (Nr+1), "mu_edges" (Nμ+1)
- root attributes for metadata (strings)
"""
function save_hist_h5(path::AbstractString, H4;
                      rmin::Real, rmax::Real, Nr::Integer,
                      mumax::Real, Nmu::Integer,
                      periodic::Bool=false, Lx=NaN, Ly=NaN, Lz=NaN)

    r_edges = collect(range(float(rmin), float(rmax); length=Int(Nr)+1))
    μ_edges = collect(range(0.0, float(mumax); length=Int(Nmu)+1))

    h5open(path, "w") do f
        f["hist"]     = H4.h
        f["r_edges"]  = r_edges
        f["mu_edges"] = μ_edges

        attrs = attributes(f)
        attrs["format"]     = "TriCo:hist4d/v1"
        attrs["axes_order"] = "r12,r23,mu12,mu13"
        attrs["rmin"]       = string(float(rmin))
        attrs["rmax"]       = string(float(rmax))
        attrs["Nr"]         = string(Int(Nr))
        attrs["mumax"]      = string(float(mumax))
        attrs["Nmu"]        = string(Int(Nmu))
        attrs["los"]        = periodic ? "z" : "true"
        attrs["periodic"]   = string(periodic)
        attrs["Lx"]         = string(float(Lx))
        attrs["Ly"]         = string(float(Ly))
        attrs["Lz"]         = string(float(Lz))
    end
    return path
end

"""
    load_hist_h5(path) -> (hist, r_edges, mu_edges, meta::Dict{String,String})

Load the histogram and metadata from an HDF5 file.
"""
function load_hist_h5(path::AbstractString)
    h5open(path, "r") do f
        hist     = read(f["hist"])
        r_edges  = read(f["r_edges"])
        mu_edges = read(f["mu_edges"])
        meta = Dict{String,String}()
        attrs = attributes(f)
        for k in keys(attrs)
            meta[string(k)] = string(attrs[k])
        end
        return hist, r_edges, mu_edges, meta
    end
end

