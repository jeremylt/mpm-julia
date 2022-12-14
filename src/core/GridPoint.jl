# ------------------------------------------------------------------------------
# Grid point
# ------------------------------------------------------------------------------

mutable struct GridPoint
    # internal state data
    isfixed::Array{Bool} # 'fixation' in different directions
    x::Array{Float64}
    m::Float64
    p::Array{Float64}
    f::Array{Float64} # total force = internal force + external force

    # constructor
    #   set all values to 0 initially
    GridPoint(x::Float64, y::Float64) =
        (new([false, false], [x, y], 0.0, zeros(2), zeros(2)))
end

# ------------------------------------------------------------------------------
