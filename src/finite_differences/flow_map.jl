using LinearAlgebra
using DifferentialEquations

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
    ∇F_fd!(dest, u, x₀, ts, δx; t₀ = 0, solver...)
Centred finite-difference approximation of the flow map gradient, vectorised on intial conditions.
"""
function ∇F_fd!(dest, u, x₀, ts, δx; t₀ = 0, solver...)
    if ndims(x₀) == 0
        # Scalar, convert to matrix
        x₀ = [x₀;;]
    elseif ndims(x₀) == 1
        # Vector, convert to matrix
        x₀ = reshape(x₀, 1, :)
    end

    n = size(x₀)[1]
    d = size(x₀)[2]
    T = maximum(ts)

    # Ensure the destination is of appropriate size
    if size(dest) != (n, length(ts))
        DimensionMismatch(
            "dest has incorrect dimensions: expected ($(n), $(length(ts))) but found $(size(dest))",
        )
    end

    # Form a star grid around each point
    # TODO: Vectorise this
    # Coordinate shifts for computing the star grid
    A = zeros(2 * d, d)
    A[1:(2 * d + 2):(2 * d^2)] .= 1
    A[2:(2 * d + 2):(2 * d^2)] .= -1

    # Index as (initial condition star point, dim)
    grid = zeros(n, 2 * d, d)
    for (i, x) in enumerate(eachrow(x₀))
        grid[i, :, :] = repeat(x', 2 * d) + δx * A
    end

    # Solve forward as an Ensemble ODEProblem
    # Advect these points forwards to the initial time
    # Mapping between grid and cartesian index i: sol[i] ⟺ grid[cld(i, 2d), mod1(i, 2d), :]
    advected = Array{Float64}(undef, n, 2 * d, d, length(ts))
    prob = ODEProblem(u, grid[1, 1, :], (t₀, T))
    # Hack the output function to directly place the solution in advected. Can ignore the output of solve.
    ensemble = EnsembleProblem(
        prob;
        prob_func = (prob, i, _) -> remake(prob; u0 = grid[cld(i, 2 * d), mod1(i, 2 * d), :]),
        output_func = function (sol, i)
            advected[cld(i, 2 * d), mod1(i, 2 * d), :, :] = sol[:, :]
            return (nothing, false)
        end,
    )
    _ = solve(ensemble, EnsembleThreads(); saveat = ts, trajectories = 2 * d * n, solver...)

    # Compute each flow map gradient, at each time step
    for i = 1:n
        for j = 1:length(ts)
            dest[i, j] = 1 / (2 * δx) * ones(d, d)
            for l = 1:d
                for m = 1:d
                    dest[i, j][l, m] *= advected[i, 2 * m - 1, l, j] - advected[i, 2 * m, l, j]
                end
            end
        end
    end
end

"""
    ∇F_fd(u, x₀, ts, δx; t₀ = 0, solver...)
Out-of-place version of ∇F_fd!
"""
function ∇F_fd(u, x₀, ts, δx; t₀ = 0, solver...)
    ∇Fs = Array{Matrix}(undef, size(x₀)[1], length(ts))
    ∇F_fd!(∇Fs, u, x₀, ts, δx; t₀ = t₀, solver...)
    return ∇Fs
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
