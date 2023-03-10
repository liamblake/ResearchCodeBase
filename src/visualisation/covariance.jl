using LinearAlgebra

using CairoMakie

"""
	bivariate_gaussian_std_dev(ax, μ, Σ; nσ = 1, colour = :black, label = nothing, kwargs...)
Plot the n standard-deviation regions of a bivariate Gaussian distribution
with mean μ and covariance matrix Σ. The number of regions plotted is specified
by nσ.
"""
function bivariate_std_dev!(ax, μ, Σ; nσ = 1, colour = :black, label = nothing, kwargs...)
    # Calculate the first two principal axes of the covariance matrix
    # These correspond to the major and minor axes of the ellipse
    evals, evecs = eigen(Σ)

    # Angle of rotation - use the principal axis
    θ = atan(evecs[2, 1], evecs[1, 1])

    # Magnitude of major and minor axes
    a, b = sqrt.(evals[1:2])

    ts = 0:0.001:(2π)

    # Plot each contour
    for n = 1:nσ
        # Parametric equations for the resulting ellipse
        # TODO: Should be a way to calculate this by operating directly on the eigenvectors
        # i.e. x = cos(θ), y = sin(θ)
        x = t -> n * (a * cos(t) * cos(θ) - b * sin(t) * sin(θ)) + μ[1]
        y = t -> n * (a * cos(t) * sin(θ) + b * sin(t) * cos(θ)) + μ[2]

        lines!(ax, x.(ts), y.(ts); color = colour, label = (n == 1) ? label : nothing, kwargs...)
    end

    # Also plot the mean
    scatter!(ax, [μ[1]], [μ[2]]; markersize = 9, color = colour, markerstrokecolor = colour)
end