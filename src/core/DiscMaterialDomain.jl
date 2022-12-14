# ------------------------------------------------------------------------------
# domain of MaterialPoints in a disc
# ------------------------------------------------------------------------------

using Printf

function creatediscmaterialdomain(
    center::Array{Float64},
    radius::Float64,
    pointsize::Float64,
    velocity::Array{Float64},
    gravity::Float64,
    density::Float64,
    E::Float64,
    ν::Float64,
)
    # adjust radius as multiple of pointsize
    radius = floor(radius / pointsize) * pointsize

    # constant parameters
    volume = pointsize * pointsize
    mass = density * volume

    # fill domain
    domain = Array{MaterialPoint}(undef, 0)
    for dy = -radius+0.5*pointsize:pointsize:radius-0.5*pointsize
        for dx = -radius+0.5*pointsize:pointsize:radius-0.5*pointsize
            if (dx^2 + dy^2 < radius^2)
                push!(
                    domain,
                    MaterialPoint(
                        mass,
                        volume,
                        [center[1] + dx, center[2] + dy],
                        velocity,
                        [0.0, gravity * mass],
                        E,
                        ν,
                    ),
                )
            end
        end
    end

    # print information
    discmass = 0.0
    for point in domain
        discmass += point.m
    end
    println("Disc domain:")
    println("  center: ", center)
    println("  radius: ", radius)
    println("  total mass: ", @sprintf("%.3f", discmass))
    println("  veloctiy (x, y): (", domain[1].v[1], ", ", domain[1].v[2], ")")
    println("  number of material points: ", length(domain))
    println("  Young's modulus: ", domain[1].E)
    println("  Poisson's ratio: ", domain[1].ν)
    println()

    # return
    domain
end

# ------------------------------------------------------------------------------
