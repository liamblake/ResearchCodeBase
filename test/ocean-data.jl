using Test

using MPhil.OceanData

@testset "Earth conversions" begin
    @test isapprox(e(WGS84), 0.0066943799901413165, atol = 1e-15)

    # Basically just make sure these functions don't error
    arc_to_meridonal(WGS84, 150.0)
    arc_to_parallel(WGS84, 150.0)
end
