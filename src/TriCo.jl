module TriCo

include("triangles.jl")
include("io.jl")
include("io_save.jl")  # defines save_hist_h5/load_hist_h5 here

export count_triangles_grid!, count_triangles_periodic_grid!,
       read_xyz_fits,
       save_hist_h5, load_hist_h5   # <- add these

end

