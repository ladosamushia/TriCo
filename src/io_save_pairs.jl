# src/io_save_pairs.jl
# NPZ save/load helpers for 2D pair histograms (Nr × Nμ).
#
# These mirror the 4D triangle helpers you already have, but operate on HistRMu
# (or any struct with a `.h :: AbstractArray{<:Integer,2}`).
#
# Usage:
#   using NPZ
#   save_pairs_npz("pairs.npz", H2)             # where H2.h is Nr×Nμ
#   save_pairs_npz("pairs_raw.npz", H2.h)       # raw 2D array
#   hist = load_pairs_npz("pairs.npz")          # returns the 2D array under "hist"
#
module IOSavePairs

using NPZ

export save_pairs_npz, load_pairs_npz

"""
    save_pairs_npz(path, H2) -> path

Save ONLY the 2D pair histogram array to a .npz file under the key `"hist"`.
`H2` is expected to have a field `.h` containing an `Nr×Nμ` integer array.
No metadata or edges are stored.
"""
function save_pairs_npz(path::AbstractString, H2)
    @assert hasfield(typeof(H2), :h) "Expected histogram object with field `.h`"
    A = getfield(H2, :h)
    @assert ndims(A) == 2 "Expected a 2D array for pairs; got $(ndims(A))D"
    NPZ.npzwrite(path, Dict("hist" => A))
    return path
end

"""
    save_pairs_npz(path, hist_array) -> path

Overload that accepts the raw 2D histogram array directly.
"""
function save_pairs_npz(path::AbstractString, hist::AbstractArray{<:Integer,2})
    NPZ.npzwrite(path, Dict("hist" => hist))
    return path
end

"""
    load_pairs_npz(path) -> hist

Load the 2D histogram array from a .npz file previously saved by `save_pairs_npz`.
Returns the array stored under the key `"hist"`.
"""
function load_pairs_npz(path::AbstractString)
    d = NPZ.npzread(path)
    haskey(d, "hist") || error("No 'hist' key found in NPZ file: $path")
    A = d["hist"]
    @assert ndims(A) == 2 "Expected a 2D array when reading pairs histogram; got $(ndims(A))D"
    return A
end

end # module

