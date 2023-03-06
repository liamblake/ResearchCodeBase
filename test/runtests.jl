using Test

@time begin
    @time include("covariance.jl")
    @time include("ocean-data.jl")
end