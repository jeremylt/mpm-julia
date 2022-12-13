# ------------------------------------------------------------------------------
# Grid point
# ------------------------------------------------------------------------------

mutable struct GridPoint
    # internal state data
    fixed::Array{Bool} # 'fixation' in different directions
    m::Float64
    x::Array{Float64}
    p::Array{Float64}
    f::Array{Float64} # total force = internal force + external force

    # constructor
    #   set all values to 0 initially
    GridPoint(x::Float64, y::Float64) =
        (new([false, false], 0.0, [x, y], zeros(2), zeros(2)))
end

# ------------------------------------------------------------------------------
