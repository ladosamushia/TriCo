using FITSIO

"""
    read_xyz_fits(path; xcol="X", ycol="Y", zcol="Z")

Open a FITS file at `path` and extract columns `xcol`, `ycol`, `zcol`
from the first binary table HDU.

Returns `(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64})`.
"""
function read_xyz_fits(path::AbstractString; xcol="X", ycol="Y", zcol="Z")
    FITS(path) do f
        hdu = f[2]  # usually HDU 2 is the binary table; adjust if needed
        X = read(hdu[xcol])
        Y = read(hdu[ycol])
        Z = read(hdu[zcol])
        return (Float64.(X), Float64.(Y), Float64.(Z))
    end
end

