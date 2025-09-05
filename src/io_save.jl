# src/io_save.jl

using NPZ

"""
    save_hist_npz(path, H4) -> path

Save ONLY the main 4D histogram array to a .npz file under the key `"hist"`.
`H4` is expected to have a field `.h` containing an `Nr×Nr×Nμ×Nμ` integer array.
No metadata, no edges are stored.
"""
function save_hist_npz(path::AbstractString, H4)
    @assert hasfield(typeof(H4), :h) "Expected histogram object with field `.h`"
    NPZ.npzwrite(path, Dict("hist" => H4.h))
    return path
end

"""
    save_hist_npz(path, hist_array) -> path

Overload that accepts the raw 4D histogram array directly.
"""
function save_hist_npz(path::AbstractString, hist::AbstractArray{<:Integer,4})
    NPZ.npzwrite(path, Dict("hist" => hist))
    return path
end

"""
    load_hist_npz(path) -> hist

Load the 4D histogram array from a .npz file previously saved by `save_hist_npz`.
Returns the array stored under the key `"hist"`.
"""
function load_hist_npz(path::AbstractString)
    d = NPZ.npzread(path)
    haskey(d, "hist") || error("No 'hist' key found in NPZ file: $path")
    return d["hist"]
end

