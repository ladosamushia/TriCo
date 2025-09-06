module TriCo

# …existing exports…
export count_triangles!, count_triangles_periodic!
export count_triangles_mixed!, count_triangles_mixed_periodic!   # ← add
export count_triangles_grid!, count_triangles_periodic_grid!

# …existing includes…
include("io.jl")
include("io_save.jl")
include("triangles.jl")
include("triangles_mixed.jl")   # ← add

end

