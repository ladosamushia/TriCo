module IOSave

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
        "los"        => periodic ? "z" : "true",
        "periodic"   => string(periodic),
        "Lx"         => string(float(Lx)),
        "Ly"         => string(float(Ly)),
        "Lz"         => string(float(Lz))
    )

    meta_json = JSON3.write(meta)                 # String
    meta_bytes = Vector{UInt8}(codeunits(meta_json))  # save as bytes

    npzwrite(path, Dict(
        "hist"       => H4.h,
        "r_edges"    => r_edges,
        "mu_edges"   => μ_edges,
        "meta_json"  => meta_bytes
    ))
    return path
end

function load_hist_npz(path::AbstractString)
    d = npzread(path)
    hist     = d["hist"]
    r_edges  = d["r_edges"]
    mu_edges = d["mu_edges"]

    meta = Dict{String,String}()
    if haskey(d, "meta_json")
        meta_str = String(d["meta_json"])
        meta_any = JSON3.read(meta_str)
        # convert to Dict{String,String}
        for (k,v) in pairs(meta_any)
            meta[string(k)] = string(v)
        end
    end
    return hist, r_edges, mu_edges, meta
end

end # module

