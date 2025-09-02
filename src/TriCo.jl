module TriCo

include("geometry.jl")
include("binning.jl")
include("grid.jl")
include("triangles.jl")
include("io.jl")
include("io_save.jl")

using .IOSave: save_hist_npz, load_hist_npz

export count_triangles!, count_triangles_periodic!,
       read_xyz_fits,
       save_hist_npz, load_hist_npz

end

