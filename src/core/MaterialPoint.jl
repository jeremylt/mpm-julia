# ------------------------------------------------------------------------------
# Material point
# ------------------------------------------------------------------------------

using LinearAlgebra

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
    σ_yield::Float64

    # stress/strain as vectors (Voigt notation)
    σ::Array{Float64}
    ε::Array{Float64}

    # plasticity
    α::Float64
    ε_plastic::Array{Float64}

    # constructor
    MaterialPoint(
        mass::Float64,
        volume::Float64,
        position::Array{Float64},
        velocity::Array{Float64},
        externalforce::Array{Float64},
        E::Float64,
        ν::Float64,
        σ_yield::Float64,
    ) = (
    # set all values to 0 initially
        new(
        mass,
        volume,
        volume,
        position,
        velocity,
        velocity * mass,
        externalforce,
        [false, false],
        zeros(2, 4),
        I(2),
        I(2),
        E,
        ν,
        σ_yield,
        zeros(3),
        zeros(3),
        0.0,
        zeros(3),
    ))
end

# ------------------------------------------------------------------------------
# Plastic strain increment
# ------------------------------------------------------------------------------

function incrementstrainplastic(materialpoint::MaterialPoint)
    dσ = zeros(3)
    dε_plastic = zeros(3)
    dα = 0.0

    μ = materialpoint.E / 2.0 / (1.0 + materialpoint.ν) # shear modulus
    λ =
        materialpoint.E * materialpoint.ν /
        ((1.0 + materialpoint.ν) * (1.0 - 2.0 * materialpoint.ν))
    κ = λ + μ

    I_2 = [1.0; 1.0; 0.0]
    I_2x2 = [1.0 1.0 0.0; 1.0 1.0 0.0; 0.0 0.0 0.0]
    I_dev = I(3) - 0.5 * I_2x2
    I_3x3 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.5] # 0.5 to convert engineering strain to physical one
    I_3x3_inv = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 2.0]

    # compute trial stress
    ε_dev = I_dev * materialpoint.ε
    s_trial = 2.0 * μ * I_3x3 * (ε_dev - materialpoint.ε_plastic)
    norm_s_trial = sqrt(s_trial[1]^2 + s_trial[2]^2 + 2 * s_trial[3]^2)
    σ_trial = κ * sum(materialpoint.ε[1] + materialpoint.ε[2]) * I_2 + s_trial

    # check yield condition
    k_1 = 0.0
    f_trial = norm_s_trial - (k_1 * materialpoint.α + materialpoint.σ_yield)

    if f_trial <= 0.0 # elastic update
        dσ = σ_trial - materialpoint.σ
    else # plastic update
        normal = s_trial / norm_s_trial
        λ = (norm_s_trial - k_1 * materialpoint.α - materialpoint.σ_yield) / (2.0 * μ + k_1)
        dα = λ
        dε_plastic = λ * I_3x3_inv * normal
        dσ =
            κ * sum(materialpoint.ε[1] + materialpoint.ε[2]) * I_2 + s_trial -
            2.0 * μ * λ * normal
        dσ -= materialpoint.σ
    end

    [dσ, dε_plastic, dα]
end

# ------------------------------------------------------------------------------
