struct SDEModel
    d::Int
    u::Function
    ∇u::Function
    σ::Function
end

struct SpatioTemporalInfo
    x₀s::AbstractVector
    ts::AbstractVector
    dt::Float64
    dx::Float64
end
