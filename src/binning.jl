# Stores cached constants for uniform linear r- and μ-binning and fast range checks.
struct TripletBinner
    rmin::Float64
    invΔr::Float64
    Nr::Int

    μmax::Float64
    invΔμ::Float64
    Nμ::Int

    rmin2::Float64
    rmax2::Float64
    μmax2::Float64
end

@inline TripletBinner(rmin::Real, rmax::Real, Nr::Integer,
                      μmax::Real, Nμ::Integer=2) = TripletBinner(
    float(rmin), Nr / (rmax - rmin), Int(Nr),
    float(μmax), Nμ / μmax, Int(Nμ),
    rmin*rmin, rmax*rmax, μmax*μmax
)

