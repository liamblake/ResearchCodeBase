module StochasticSensitivity
include("structs.jl")
include("covariance.jl")

export SDEModel, SpatioTemporalInfo
export compute_Σ!, compute_Σ, star_grid, ∇F

end