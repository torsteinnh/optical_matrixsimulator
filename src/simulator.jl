module simulator

export matrixcore, fresnelltools, analyticalmaterials, cmtgratings, utilities, materials

include("utilities.jl")
include("analyticalmaterials.jl")
include("materials.jl")
include("matrixcore.jl")
include("fresnelltools.jl")
# include("cmtgratings.jl")

end # simulator
