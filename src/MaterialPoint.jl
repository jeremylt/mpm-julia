# ------------------------------------------------------------------------------
# Material point
# ------------------------------------------------------------------------------

mutable struct MaterialPoint
    # internal state data
    m::Float64
    V_0::Float64
    V::Float64
    x::Array{Float64} # centroid
    v::Array{Float64}
    p::Array{Float64} # momentum
    f_external::Array{Float64}
    restraint::Array{Bool}
    corner::Array{Float64}

    # deformation gradient
    F::Array{Float64}
    dF::Array{Float64}

    # material parameters
    E::Float64
    ν::Float64

    # stress/strain as vectors (Voigt notation)
    σ::Array{Float64}
    ε::Array{Float64}

    # constructor
    MaterialPoint(
        mass::Float64,
        volume::Float64,
        position::Array{Float64},
        velocity::Array{Float64},
        externalforce::Array{Float64},
        E::Float64,
        ν::Float64,
    ) = (
    # set all values to 0 initially
        new(
        mass,
        volume,
        volume,
        position,
        velocity,
        zeros(2),
        externalforce,
        [false, false],
        zeros(2, 4),
        eye(2, 2),
        eye(2, 2),
        E,
        ν,
        zeros(3),
        zeros(3),
    ))
end

# ------------------------------------------------------------------------------
