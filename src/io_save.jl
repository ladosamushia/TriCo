module IOSave

using NPZ

"""
    save_hist_npz(path, H4; rmin, rmax, Nr, mumax, Nmu, periodic=false, Lx=NaN, Ly=NaN, Lz=NaN)

Save histogram `H4.h` to a `.npz` file with:
- dataset: `hist`      (4D array of counts, shape Nr×Nr×Nμ×Nμ)
- dataset: `r_edges`   (bin edges for r)
- dataset: `mu_edges`  (bin edges for μ)
- dataset: `meta_keys`, `meta_values` (simple metadata dict)
"""
function save_hist_npz(path::AbstractString, H4;
                       rmin::Real, rmax::Real, Nr::Integer,
                       mumax::Real, Nmu::Integer,
                       periodic::Bool=false, Lx=NaN, Ly=NaN, Lz=NaN)

    r_edges = collect(range(float(rmin), float(rmax); length=Int(Nr)+1))
    μ_edges = collect(range(0.0, float(mumax); length=Int(Nmu)+1))

    meta = Dict(
        "format"     => "TriCo:hist4d/v1",
        "axes_order" => "r12,r23,mu12,mu13",
        "rmin"       => string(float(rmin)),
        "rmax"       => string(float(rmax)),
        "Nr"         => string(Int(Nr)),
        "mumax"      => string(float(mumax)),
        "Nmu"        => string(Int(Nmu)),
        "los"        => "z",
        "periodic"   => string(periodic),
        "Lx"         => string(float(Lx)),
        "Ly"         => string(float(Ly)),
        "Lz"         => string(float(Lz))
    )

    npzwrite(path, Dict(
        "hist"        => H4.h,
        "r_edges"     => r_edges,
        "mu_edges"    => μ_edges,
        "meta_keys"   => collect(keys(meta)),
        "meta_values" => collect(values(meta))
    ))
    return path
end

"""
    load_hist_npz(path) -> (hist, r_edges, mu_edges, meta::Dict)

Load histogram and metadata back from `.npz`.
"""
function load_hist_npz(path::AbstractString)
    d = npzread(path)
    hist = d["hist"]
    r_edges = d["r_edges"]
    mu_edges = d["mu_edges"]
    meta = Dict(String.(d["meta_keys"]) .=> String.(d["meta_values"]))
    return hist, r_edges, mu_edges, meta
end

end # module

