module TriCo

# …existing exports…
export count_triangles!, count_triangles_periodic!
export count_triangles_mixed!, count_triangles_mixed_periodic!   # ← add
export count_triangles_grid!, count_triangles_periodic_grid!
export count_pairs_grid!
export count_pairs_cross_grid!
export save_pairs_npz, load_pairs_npz   # if you want them public

# …existing includes…
include("io.jl")
include("io_save.jl")
include("pairs_utils.jl")
include("triangles.jl")
include("triangles_mixed.jl")   # ← add
include("pairs.jl")
include("pairs_mixed.jl")
include("io_save_pairs.jl")
end

