using Test

@time begin
    @time include("numerics.jl")
    @time include("covariance.jl")
    @time include("ocean-data.jl")
end