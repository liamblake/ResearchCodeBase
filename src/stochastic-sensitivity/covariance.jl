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
    dt::Real,
    dx::Real;
    method::String = "fd",
    ∇u::Function = nothing,
    ode_solver = Euler(),
)
    d = length(x₀)

    if !(method in ["fd", "ode", "eov"])
        ArgumentError("method must be one of \"fd\", \"ode\" or \"eov\", got $(method).")
    end

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
    else
        # Calculate the flow map gradients by solving the equation of variations directly
        if method == "eov"
            if ∇u === nothing
                ArgumentError("Must specify ∇u to use the Equation of Variations (eov) method.")
            end

            ∇u_F = t -> ∇u(det_sol(t), t)
            ∇Fs = Vector{Matrix{Float64}}(undef, length(ts))
            ∇F_eov!(∇Fs, ∇u_F, d, t₀, T, dt)
        elseif method == "fd"
            # Use the star grid method to approximate the flow map

            # Form the star grid around the initial position
            star = star_grid(x₀, dx)

            # Advect these points forwards to the initial time
            prob = ODEProblem(velocity, star[1, :], (t₀, T))
            ensemble =
                EnsembleProblem(prob; prob_func = (prob, i, _) -> remake(prob; u0 = star[i, :]))
            sol = solve(
                ensemble,
                ode_solver,
                EnsembleThreads();
                dt = dt,
                dtmax = dt,
                saveat = ts,
                trajectories = 2 * d,
            )

            # Permute the dimensions of the ensemble solution so star_values is indexed
            # as (timestep, gridpoint, coordinate).
            star_values = Array{Float64}(undef, length(sol[1]), 2 * d, d)
            permutedims!(star_values, Array(sol), [2, 3, 1])

            # Approximate the flow map gradient at each time step
            ∇Fs = ∇F.(eachslice(star_values; dims = 1), d, dx)
        end

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
    end

    return w, Σ
end