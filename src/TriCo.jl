module TriCo

export count_triangles!, HistR12R23Mu12Mu13
export count_triangles_periodic!, HistR12R23Mu12Mu13


include("geometry.jl")   # dist² and μ² kernels
include("binning.jl")    # binner config (cached constants)
include("grid.jl")       # cell-linked list grid
include("triangles.jl")  # main hot loop

end

