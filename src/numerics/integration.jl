"""
iterative_trapezoidal!(dest::AbstractVector, ts::AbstractVector, integrand::AbstractVector)

Approximate a definite integral of the form ∫f(τ)dτ from τ = 0 to τ = t. using the trapezoidal
rule. The integral is evaluated iteratively, saving each step. The integrand input contains
evaluations f(tᵢ) corresponding to each tᵢ in ts. The computation is stored in dest, which should
have the same length as both ts and integrand.
"""
function iterative_trapezoidal!(dest::AbstractVector, ts::AbstractVector, integrand::AbstractVector)
    n = length(integrand)
    if length(ts) != n
        DimensionMismatch(
            "ts and integrand evaluations must be of same size, got $(length(ts)) and $(n).",
        )
    end

    if length(dest) != n
        DimensionMismatch(
            "Destination must match size of input, expected $(n) but got $(length(dest))",
        )
    end

    # Compute the trapezoidal area, saving at every step
    dest[1] = zero(integrand[1])
    @inbounds for i = 2:n
        dest[i] = dest[i - 1] + 0.5 * abs(ts[i] - ts[i - 1]) * (integrand[i - 1] + integrand[i])
    end
end