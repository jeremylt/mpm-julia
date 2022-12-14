# ------------------------------------------------------------------------------
# domain of MaterialPoints in a box
# ------------------------------------------------------------------------------

using Printf

function createboxmaterialdomain(
    lowerleftcorner::Array{Float64},
    width::Float64,
    height::Float64,
    pointsize::Float64,
    velocity::Array{Float64},
    gravity::Float64,
    density::Float64,
    E::Float64,
    ν::Float64,
    σ_yield::Float64,
)
    # adjust width and height as multiple of pointsize
    width = floor(width / pointsize) * pointsize
    height = floor(height / pointsize) * pointsize

    # constant parameters
    volume = pointsize * pointsize
    mass = density * volume

    # fill domain
    domain = Array{MaterialPoint}(undef, 0)
    for dy = 0.0:pointsize:height
        for dx = 0.0:pointsize:width
            push!(
                domain,
                MaterialPoint(
                    mass,
                    volume,
                    [lowerleftcorner[1] + dx, lowerleftcorner[2] + dy],
                    velocity,
                    [0.0, gravity * mass],
                    E,
                    ν,
                    σ_yield,
                ),
            )
        end
    end

    # print information
    boxmass = 0.0
    for point in domain
        boxmass += point.m
    end
    println("Box domain:")
    println("  lower left corner: ", lowerleftcorner)
    println("  width (x): ", width)
    println("  height (y): ", height)
    println("  total mass: ", @sprintf("%.3f", boxmass))
    println("  veloctiy (x, y): (", domain[1].v[1], ", ", domain[1].v[2], ")")
    println("  number of material points: ", length(domain))
    println("  Young's modulus: ", domain[1].E)
    println("  Poisson's ratio: ", domain[1].ν)
    println("  Yield stress: ", domain[1].σ_yield)
    println()

    # return
    domain
end

# ------------------------------------------------------------------------------
