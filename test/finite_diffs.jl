using Test

using MPhil.FiniteDiffs

@testset "∇F linear" begin
    """
    Tests the calculation of ∇F with a linear system.
    """
    # Define velocity
    A = [1.0 0.0; 0.0 2.3]
    u = (x, _, _) -> A * x
    ∇u = (_, _) -> A

    xs = [1.1 -0.2; 5.2 0.1; 0.5 -9.0]
    ts = [0.2, 0.5, 0.9]

    # Exact forms of terms involved
    Fe = (x, t) -> [exp(A[1, 1] * t) * x[1]; exp(A[2, 2] * t) * x[2]]
    ∇Fe = (x, t) -> [exp(A[1, 1] * t) 0.0; 0.0 exp(A[2, 2] * t)]

    ∇Fes = Array{Matrix}(undef, size(xs)[1], length(ts))
    for (i, x) in enumerate(eachrow(xs))
        for (j, t) in enumerate(ts)
            ∇Fes[i, j] = ∇Fe(x, t)
        end
    end

    for method in ["fd"]#, "eov"]
        ∇Fs = ∇F(u, xs, ts; δx = 0.001, method = method, ∇u = ∇u)
        @test isapprox(∇Fs, ∇Fes, atol = 1e-3)
    end
end
