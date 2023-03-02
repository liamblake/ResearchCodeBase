using Test

using MPhil.StochasticSensitivity

@testset "Σ_calculation errors" begin end

@testset "OU calculations" begin
    """
    Tests the calculation of Σ with an OU process, for which Σ can be computed exactly. See the
    supplementary materials for more details.
    """
    # Define OU process
    A = [1.0 0.0; 0.0 2.3]
    σ = [1.0 0.5; -0.2 1.4]
    u = (x, _, _) -> A * x
    ∇u = (x, t) -> A
    σf = (_, _, _) -> σ

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
    for method in ["∇F", "ode"]
        # Test full w, Σ calculation
        w, Σ = Σ_calculation(u, σf, x, 0, t, 0.001; method = method, δx = 0.001, ∇u = ∇u)

        # TODO: Fix tolerances
        @test isapprox(w, Fe, atol = 1e-0)
        @test_skip isapprox(Σ, Σe, atol = 1e-0)
    end
end
