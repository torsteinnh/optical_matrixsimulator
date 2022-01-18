module simulator

export matrixcore, fresnelltools, analyticalmaterials, cmtgratings, utilities, materials

include("utilities.jl")
include("analyticalmaterials.jl")
include("ema.jl")
include("materials.jl")
include("matrixcore.jl")
include("fresnelltools.jl")

end # simulator
