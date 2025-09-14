module TriCo

# ================= Exports =================
export count_triangles!, count_triangles_periodic!
export count_triangles_mixed!, count_triangles_mixed_periodic!
export count_triangles_grid!, count_triangles_periodic_grid!
export count_pairs_grid!, count_pairs_cross_grid!
export save_pairs_npz, load_pairs_npz

# ================= Includes =================
include("io.jl")
include("io_save.jl")
include("pairs_utils.jl")
include("triangles.jl")
include("triangles_mixed.jl")
include("pairs.jl")
include("pairs_mixed.jl")
include("io_save_pairs.jl")

# ================= GPU Support (optional) =================
try
    include("triangles_gpu.jl")
    using .TrianglesGPU
    export TrianglesGPU
    # If you want top-level exports as well:
    # export count_triangles_grid_gpu!, count_triangles_periodic_grid_gpu!
catch e
    @warn "GPU support not loaded. CUDA may be unavailable or CUDA.jl not installed." exception=(e, catch_backtrace())
end

end # module TriCo

