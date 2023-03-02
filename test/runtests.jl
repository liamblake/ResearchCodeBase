using Test

@time begin
    @time include("finite_diffs.jl")
    @time include("covariance.jl")
end