using Test

using DifferentialEquations

using MPhil.StochasticSensitivity

@testset "star grid" begin
    """
    Check that the star grid computation is correct, given a point x and perturbation δx.
    """
    x = [1.0, 2.0, 3.0]
    δx = 0.01

    sg = star_grid(x, δx)

    @test sg == [1.01 2.0 3.0; 0.99 2.0 3.0; 1.0 2.01 3.0; 1.0 1.99 3.0; 1.0 2.0 3.01; 1.0 2.0 2.99]
end

@testset "gradient approximation" begin
    """
    Check that the finite-difference approximation of ∇F gives the right calculation, given values of
    the stargrid.
    """
    δx = 1.0
    star = [[2.0 0.0]; [0.8 0.2]; [1.6 0.6]; [1.4 -0.2]]

    expected = [[0.6 0.1]; [-0.1 0.4]]

    @test isapprox(∇F(star, 2, δx), expected, atol = 1e-10)
end

@testset "Σ_calculation errors" begin end

@testset "OU calculations" begin
    """
    Tests the calculation of Σ with an OU process, for which Σ can be computed exactly. See the
    supplementary materials for more details.
    """
    # Define OU process
    A = [1.0 0.0; 0.0 2.3]
    σ = [1.0 0.5; -0.2 1.4]
    u = (x, _) -> A * x
    ∇u = (x, t) -> A
    σf = (_, _) -> σ

    x = [1.1, -0.2]
    t = 1.5

    # Exact forms of terms involved
    Fe = [exp(A[1, 1] * t) * x[1]; exp(A[2, 2] * t) * x[2]]
    ∇Fe = [exp(A[1, 1] * t) 0.0; 0.0 exp(A[2, 2] * t)]
    Σe = [
        (σ[1, 1]^2 + σ[1, 2]^2) / (2 * A[1, 1])*(exp(2 * A[1, 1] * t) - 1) (σ[1, 1] * σ[2, 1] + σ[1, 2] * σ[2, 2]) / (A[1, 1] + A[2, 2])*(exp((A[1, 1] + A[2, 2]) * t) - 1)
        (σ[1, 1] * σ[2, 1] + σ[1, 2] * σ[2, 2]) / (A[1, 1] + A[2, 2])*(exp((A[1, 1] + A[2, 2]) * t) - 1) (σ[2, 1]^2 + σ[2, 2]^2) / (2 * A[2, 2])*(exp(2 * A[2, 2] * t) - 1)
    ]

    # Test each method
    for method in ["flow map", "backwards flow map", "ode rk4"]
        # Test full w, Σ calculation
        w, Σ, _, _ = compute_Σ(
            SDEModel(2, u, ∇u, σf),
            SpatioTemporalInfo([x], 0:0.001:t, 0.001, 0.001),
            method,
        )

        @test w[:, 1] == [x]
        @test isapprox(w[1, end], Fe, atol = 1e-1)
        @test Σ[:, 1] == [zeros(2, 2)]
        @test isapprox(Σ[1, end], Σe, atol = 1e-1) broken = (method == "backwards flow map")
    end
end

@testset "GBM 1D" begin
    # Define GBM process
    a = 1.1
    b = 0.29
    u = (x, _) -> a * x
    ∇u = (x, t) -> [a]
    σ = (x, _) -> b * x

    x = 1.1
    t = 1.5

    # Exact forms of terms involved
    Fe = [x * exp(a * t)]
    ∇Fe = [exp(a * t);;]
    Σe = [t * b^2 * x^2 * exp(2 * a * t);;]

    # Test each method
    for method in ["flow map", "backwards flow map", "ode rk4"]
        # Test full w, Σ calculation
        w, Σ, _, _ = compute_Σ(
            SDEModel(1, u, ∇u, σ),
            SpatioTemporalInfo([[x]], 0:0.001:t, 0.001, 0.001),
            method,
        )

        @test w[:, 1] == [[x]]
        @test isapprox(w[1, end], Fe, atol = 1e-2)
        @test Σ[:, 1] == [[0.0;;]]
        @test isapprox(Σ[1, end], Σe, atol = 1e-2)
    end
end