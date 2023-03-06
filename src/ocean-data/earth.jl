"""
Describes an ellipsoidal model of the Earth, with semi-major axis a and inverse flattening inv_f.
"""
struct EarthModel
    a::Float64
    inv_f::Float64
end

"""
    e(self::EarthModel)::Float64

Calculate the eccentricity of a symmetric cross-sectional ellipse of an EarthModel.
"""
function e(self::EarthModel)::Float64
    f = 1.0 / self.inv_f
    return 2 * f - f^2
end

WGS84 = EarthModel(6378187.0, 298.257223563);

"""
    arc_to_meridonal(self::EarthModel, φ)

Compute the length (in degrees) of a 1 metre change in latitude, at latitude φ.
"""
function arc_to_meridonal(self::EarthModel, φ)
    e² = e(self)^2
    return @. 180 / pi * (1 - e²) / self.a * (1 - e² * sind(φ)^2)^(3 / 2)
end
"""
    arc_to_meridonal(self::EarthModel, φ)

Compute the length (in degrees) of a 1 metre change in longitude, at latitude φ.
"""
function arc_to_parallel(self::EarthModel, φ)
    e² = e(self)^2
    return @. 180 / pi * (1 - e²) / self.a * (1 - e² * sind.(φ)^2)^(1 / 2) / cosd.(φ)
end