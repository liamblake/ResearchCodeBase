using Test

using MPhil.Numerics

@testset "trapezoidal rule" begin
    ts = [0, 0.01, 0.02]
    dest = Vector{Float64}(undef, length(ts))
    integrand = @. ts * cos(ts^2)

    iterative_trapezoidal!(dest, ts, integrand)
    truth = t -> 0.5 * sin(t^2)

    @test isapprox(dest, truth.(ts), atol = 1e-6)

    # Should also work with a matrix state
    dest2 = Vector{Matrix}(undef, length(ts))
    f = t -> [t*cos(t^2) exp(t); 0.0 1.0]
    integrand2 = f.(ts)

    truth2 = t -> [truth(t) exp(t)-1; 0.0 t]
    iterative_trapezoidal!(dest2, ts, integrand2)

    @test isapprox(dest2, truth2.(ts), atol = 1e-6)
end
