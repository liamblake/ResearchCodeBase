using LinearAlgebra

using DifferentialEquations
using Parameters

"""
compute_Σ!(
    w_dest,
    Σ_dest,
    eval_dest,
    θ_dest,
    model::SDEModel,
    st_info::SpatioTemporalInfo,
    method;
    kwargs...
)

Calculate a deterministc trajectory and corresponding covariance matrices, at times specified by
ts, using one of several different methods.

Available methods are:
    - flow map: Use the forwards version of the flow-map integral expression for Σ.
    - backwards flow map: Use the backwards version of the flow map integral expression for Σ.
    - ode rk4: Solve the ODE that Σ satisfies directly with the order-4 Runge-Kutta scheme.
"""
function compute_Σ!(
    w_dest,
    Σ_dest,
    eval_dest,
    θ_dest,
    model::SDEModel,
    st_info::SpatioTemporalInfo,
    method;
    kwargs...,
)
    if method == "flow map"
        Σ_flow_map!(
            w_dest,
            Σ_dest,
            eval_dest,
            θ_dest,
            model,
            st_info;
            backwards = false,
            eov = false,
            kwargs...,
        )
    elseif method == "backwards flow map"
        Σ_flow_map!(
            w_dest,
            Σ_dest,
            eval_dest,
            θ_dest,
            model,
            st_info;
            backwards = true,
            eov = false,
            kwargs...,
        )
    elseif method == "flow map eov"
        Σ_flow_map!(
            w_dest,
            Σ_dest,
            eval_dest,
            θ_dest,
            model,
            st_info;
            backwards = true,
            eov = true,
            kwargs...,
        )
    elseif method == "ode rk4"
        Σ_ode_rk4!(w_dest, Σ_dest, eval_dest, θ_dest, model, st_info; kwargs...)
    else
        ArgumentError("Invalid method for computing Σ... got $(method)")
    end

    nothing
end

"""

Out-of-place version of compute_Σ!, allocating and returning results.
"""
function compute_Σ(model::SDEModel, st_info::SpatioTemporalInfo, method; kwargs...)
    w_dest = Array{Vector}(undef, length(st_info.x₀s), length(st_info.ts))
    Σ_dest = Array{Matrix}(undef, length(st_info.x₀s), length(st_info.ts))
    eval_dest = Array{Vector}(undef, length(st_info.x₀s), length(st_info.ts))
    θ_dest = Array{Vector}(undef, length(st_info.x₀s), length(st_info.ts))

    compute_Σ!(w_dest, Σ_dest, eval_dest, θ_dest, model, st_info, method; kwargs...)

    return w_dest, Σ_dest, eval_dest, θ_dest
end

"""
	star_grid(x, δx)
Construct a star grid around a point x, for calculating a finite difference
approximation of a spatial derivative. The number of points in the grid is
determined from the dimension of x; if the length of x is n, then 2n points are
calculated.
"""
function star_grid(x::AbstractVector, δx::Real)
    n = length(x)

    # Matrix of coordinate shifts
    # ⌈1  0  0  0  ⋯ 0⌉
    # |-1 0  0  0  ⋯ 0|
    # |0  1  0  0  ⋯ 0|
    # |0  -1 0  0  ⋯ 0|
    # |⋮     ⋱       ⋮|
    # |0  ⋯     0  ⋯ 1|
    # ⌊0  ⋯     0  ⋯ 0⌋
    A = zeros(2 * n, n)
    A[1:(2 * n + 2):(2 * n^2)] .= 1
    A[2:(2 * n + 2):(2 * n^2)] .= -1

    return repeat(x', 2 * n) + δx * A
end

"""
    ∇F(star_values, n, δx)
Approximate the flow map gradient with a centered finite-difference approximation, given a star grid of values.
"""
function ∇F(star_values, n, δx)
    return 1 / (2 * δx) * (star_values[1:2:(2 * n), :] - star_values[2:2:(2 * n), :])'
end

"""
    ∇F_eov!(dest, ∇u, d, t₀, T, dt)
Calculate the flow map gradient by directly solving the equation of variations,
given the corresponding trajectory. The equation of variations is
	∂∇F/∂t = ∇u(F(t), t)∇F.
The gradient of the flow map is taken with respect to the initial condition at time
t₀. A vector of matrices, corresponding to the flow map gradient at time steps at dt,
is placed in the preallocated dest vector.
"""
function ∇F_eov!(dest, ∇u, d, t₀, T, dt)
    # Inplace definition of the ODE
    function rate(x, _, t)
        dx = ∇u(t) * x
        return dx
    end

    Id = zeros(d, d)
    Id[diagind(Id)] .= 1.0
    u₀ = Id

    prob = ODEProblem(rate, u₀, (t₀, T))
    dest[:] = solve(prob; saveat = dt).u
end

"""
	Σ_calculation(model, x₀, t₀, T, dt)
Calculate the deviation covariance matrix Σ with an in-place specification of the velocity field.
"""
function Σ_flow_map!(
    w_dest,
    Σ_dest,
    eval_dest,
    θ_dest,
    model::SDEModel,
    st_info::SpatioTemporalInfo;
    ode_solver = RK4(),
    backwards = false,
    eov = false,
    solver_opts...,
)
    @unpack d, u, σ = model
    @unpack x₀s, ts, dt, dx = st_info
    n = size(x₀s)[1]

    t₀ = minimum(ts)
    T = maximum(ts)

    @inbounds for (j, x₀) in enumerate(x₀s)
        # Generate the required flow map data
        # First, advect the initial condition forward to obtain the final position
        prob = ODEProblem((x, _, t) -> u(x, t), x₀, (t₀, T))
        det_sol = solve(prob, ode_solver; dt = dt, saveat = ts, adaptive = false, solver_opts...)
        w_dest[j, :] = det_sol.u

        if backwards
            # Form the star grid around the final position
            star = star_grid(w_dest[j, end], dx)

            # Advect these points forwards to the initial time
            prob = ODEProblem((x, _, t) -> u(x, t), star[1, :], (T, t₀))
            ensemble =
                EnsembleProblem(prob; prob_func = (prob, i, _) -> remake(prob; u0 = star[i, :]))
            sol = solve(
                ensemble,
                ode_solver,
                EnsembleThreads();
                trajectories = 2 * d,
                dt = dt,
                saveat = reverse(ts),
                adaptive = false,
                solver_opts...,
            )
        else
            # Form the star grid around the initial position
            star = star_grid(x₀, dx)

            # Advect these points forwards to the initial time
            prob = ODEProblem((x, _, t) -> u(x, t), star[1, :], (t₀, T))
            ensemble =
                EnsembleProblem(prob; prob_func = (prob, i, _) -> remake(prob; u0 = star[i, :]))
            sol = solve(
                ensemble,
                ode_solver,
                EnsembleThreads();
                trajectories = 2 * d,
                dt = dt,
                saveat = ts,
                adaptive = false,
                solver_opts...,
            )
        end

        # Permute the dimensions of the ensemble solution so star_values is indexed
        # as (timestep, gridpoint, coordinate).
        star_values = Array{Float64}(undef, length(sol[1]), 2 * d, d)
        permutedims!(star_values, Array(sol), [2, 3, 1])

        # Approximate the flow map gradient at each time step
        ∇Fs = ∇F.(eachslice(star_values; dims = 1), d, dx)
        if backwards
            ∇Fs[end] = Matrix{Float64}(I, d, d)
            reverse!(∇Fs)
        else
            ∇Fs[1] = Matrix{Float64}(I, d, d)
        end

        # Evaluate σ along the trajectory
        σs = σ.(det_sol.(ts), ts)

        # Evaluate the integrand at each time point.
        # See Theorem 2.2 for the integral form of Σ.
        Ks = inv.(∇Fs) .* σs

        integrand = Ks .* transpose.(Ks)

        # Approximate the inner integral using the trapezoidal rule
        Σ_dest[j, 1] = zeros(d, d)
        prev_val = Σ_dest[j, 1]
        eval_dest[j, 1] = zeros(d)
        θ_dest[j, 1] = zeros(d)
        for i = 2:length(ts)
            prev_val += 0.5 * abs(ts[i] - ts[i - 1]) * (integrand[i - 1] + integrand[i])
            Σ_dest[j, i] = prev_val

            if !backwards
                Σ_dest[j, i] = ∇Fs[i] * Σ_dest[j, i] * transpose(∇Fs[i])
            end

            # Compute stochastic sensitivity
            E = eigen(Σ_dest[j, i]; sortby = λ -> (-real(λ), -imag(λ)))
            eval_dest[j, i] = E.values
            if d == 2
                θ_dest[j, i] = atan.(E.vectors[2, :], E.vectors[1, :])
            end
        end
    end
    nothing
end

function Σ_ode_rk4!(w_dest, Σ_dest, eval_dest, θ_dest, model::SDEModel, st_info::SpatioTemporalInfo)
    @unpack d, u, ∇u, σ = model
    @unpack ts, x₀s, dt = st_info

    # Preallocate temp values
    kw1 = Vector{Float64}(undef, d)
    kw2 = Vector{Float64}(undef, d)
    kw3 = Vector{Float64}(undef, d)
    kw4 = Vector{Float64}(undef, d)
    kΣ1 = Matrix{Float64}(undef, d, d)
    kΣ2 = Matrix{Float64}(undef, d, d)
    kΣ3 = Matrix{Float64}(undef, d, d)
    kΣ4 = Matrix{Float64}(undef, d, d)

    w_dest[:, 1] = x₀s
    @inbounds for i = 1:(size(Σ_dest)[1])
        Σ_dest[i] = zeros(d, d)
        eval_dest[i, 1] = zeros(d)
        θ_dest[i, 1] = zeros(d)

        kw1 .= 0.0
        kw2 .= 0.0
        kw3 .= 0.0
        kw4 .= 0.0
        kΣ1 .= 0.0
        kΣ2 .= 0.0
        kΣ3 .= 0.0
        kΣ4 .= 0.0

        @inbounds for j = 2:length(ts)
            # RK4 method
            kw1 .= u(w_dest[i, j - 1], ts[j - 1])
            kw2 .= u(w_dest[i, j - 1] + dt / 2 * kw1, ts[j - 1] + dt / 2)
            kw3 .= u(w_dest[i, j - 1] + dt / 2 * kw2, ts[j - 1] + dt / 2)
            kw4 .= u(w_dest[i, j - 1] + dt * kw3, ts[j - 1] + dt)
            w_dest[i, j] = w_dest[i, j - 1] + dt / 6 * (kw1 + 2 * kw2 + 2 * kw3 + kw4)

            kΣ1 .=
                ∇u(w_dest[i, j - 1], ts[j - 1]) * Σ_dest[i, j - 1] +
                Σ_dest[i, j - 1] * transpose(∇u(w_dest[i, j - 1], ts[j - 1])) +
                σ(w_dest[i, j - 1], ts[j - 1]) * transpose(σ(w_dest[i, j - 1], ts[j - 1]))

            kΣ2 .=
                ∇u(w_dest[i, j - 1] + dt / 2 * kw1, ts[j - 1] + dt / 2) *
                (Σ_dest[i, j - 1] + dt / 2 * kΣ1) +
                (Σ_dest[i, j - 1] + dt / 2 * kΣ1) *
                transpose(∇u(w_dest[i, j - 1] + dt / 2 * kw1, ts[j - 1] + dt / 2)) +
                σ(w_dest[i, j - 1] + dt / 2 * kw1, ts[j - 1] + dt / 2) *
                transpose(σ(w_dest[i, j - 1] + dt / 2 * kw1, ts[j - 1] + dt / 2))

            kΣ3 .=
                ∇u(w_dest[i, j - 1] + dt / 2 * kw2, ts[j - 1] + dt / 2) *
                (Σ_dest[i, j - 1] + dt / 2 * kΣ2) +
                (Σ_dest[i, j - 1] + dt / 2 * kΣ2) *
                transpose(∇u(w_dest[i, j - 1] + dt / 2 * kw2, ts[j - 1] + dt / 2)) +
                σ(w_dest[i, j - 1] + dt / 2 * kw2, ts[j - 1] + dt / 2) *
                transpose(σ(w_dest[i, j - 1] + dt / 2 * kw2, ts[j - 1] + dt / 2))

            kΣ4 .=
                ∇u(w_dest[i, j - 1] + dt * kw3, ts[j - 1] + dt) * (Σ_dest[i, j - 1] + dt * kΣ3) +
                (Σ_dest[i, j - 1] + dt * kΣ3) *
                transpose(∇u(w_dest[i, j - 1] + dt * kw3, ts[j - 1] + dt)) +
                σ(w_dest[i, j - 1] + dt * kw3, ts[j - 1] + dt) *
                transpose(σ(w_dest[i, j - 1] + dt * kw3, ts[j - 1] + dt))

            Σ_dest[i, j] = Σ_dest[i, j - 1] + dt / 6 * (kΣ1 + 2 * kΣ2 + 2 * kΣ3 + kΣ4)

            # Compute stochastic sensitivity measures
            E = eigen(Σ_dest[i, j]; sortby = λ -> (-real(λ), -imag(λ)))
            eval_dest[i, j] = E.values
            if d == 2
                θ_dest[i, j] = atan.(E.vectors[2, :], E.vectors[1, :])
            end
        end
    end
    nothing
end
