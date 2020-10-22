module AnimatedOptimization

using Distributions, Plots, ForwardDiff, LinearAlgebra,
  Convex, ECOS

export minrandomsearch,
  graddescent,
  newton,
  interiorpoint,
  sequentialquadratic,
  slqp
  

include("heuristic_optimizers.jl")
include("smooth_optimizers.jl")
include("constrained_optimizers.jl")


end # module AnimatedOptimization
