using LinearAlgebra

using DifferentialEquations

using ..FiniteDiffs

"""
	Σ_calculation(model, x₀, t₀, T, dt)
Calculate the deviation covariance matrix Σ with an in-place specification of the velocity field.
"""
function Σ_calculation(
    velocity::Function,
    σ::Function,
    x₀::AbstractVector,
    t₀::Real,
    T::Real,
    dt::Real;
    method::String = "∇F",
    ∇F_method::String = "fd",
    δx::Real = nothing,
    ∇u::Function = nothing,
    ode_solver = Euler(),
    ∇F_kwargs...,
)
    # TODO: Vectorise across the initial condition
    d = length(x₀)

    ts = t₀:dt:T
    if last(ts) < T
        ts = range(t₀; stop = T, length = length(ts) + 1)
    end

    # Generate the required flow map data
    # First, advect the initial condition forward to obtain the final position
    prob = ODEProblem(velocity, x₀, (t₀, T))
    det_sol = solve(prob, ode_solver; dt = dt, dtmax = dt, saveat = ts)
    w = last(det_sol)

    if method == "ode"
        # Calculate the covariance by directly solving the differential equation
        if ∇u === nothing
            ArgumentError("Must specify ∇u to use the ODE method.")
        end

        # Solve for the flow map and the covariance matrix as a joint system
        function joint_F_Σ(x, _, t)
            traj_∇u = ∇u(x[1], t)
            traj_σ = σ(x[1], nothing, t)

            du = [velocity(x[1], nothing, t), traj_∇u * x[2] + x[2] * traj_∇u' + traj_σ * traj_σ']
            return du
        end

        prob = ODEProblem(joint_F_Σ, [x₀, zeros(d, d)], (t₀, T))
        sol = solve(prob, ode_solver; dt = dt, dtmax = dt, save_everystep = false)

        Σ = sol[2][2]

    elseif method == "∇F"
        if ∇F_method == "fd" && δx === nothing
            ArgumentError("Must specify δx to use ∇F method with finite-difference.")
        end

        ∇Fs = Array{Float64}(undef, 1, length(ts))
        ∇F!(∇Fs, velocity, x₀, ts; δx = δx, t₀ = 0, method = ∇F_method, ∇u = ∇u, ∇F_kwargs...)

        # Evaluate σ along the trajectory
        σs = σ.(det_sol.(ts), nothing, ts)

        # Evaluate the integrand at each time point.
        # See Theorem 2.2 for the integral form of Σ.
        Ks = [last(∇Fs)] .* inv.(∇Fs) .* σs
        integrand = Ks .* transpose.(Ks)

        # Approximate an integral given discrete evaluations, using the composite Simpson's rule
        # Assume the data is equidistant with respect to the variable being integrate.
        feven = @view integrand[2:2:end]
        fodd = @view integrand[1:2:end]

        Σ = dt / 3 * (integrand[1] + 2 * sum(feven) + 4 * sum(fodd) + last(integrand))

    else
        ArgumentError("method must be one of \"∇F\" or \"ode\" but got $(method).")
    end

    return w, Σ
end