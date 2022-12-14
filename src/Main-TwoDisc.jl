# ------------------------------------------------------------------------------
# 2D example - plate and disc
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# object definitions
# ------------------------------------------------------------------------------

include("MaterialPoint.jl")
include("DiscMaterialDomain.jl")
include("GridPoint.jl")
include("Grid.jl")
include("Basis.jl")

# ------------------------------------------------------------------------------
# driver
# ------------------------------------------------------------------------------

function main()
    # physical constants
    g = 0.0
    ρ = 1000.0 # density
    E = 1000.0 # Young's modulus
    ν = 0.3    # Poisson's ratio

    # setup grid
    println("--------- Background grid ---------")
    numcells = 20
    grid = Grid(1.0, 1.0, numcells + 1, numcells + 1)

    # setup disc parameters
    radius = 0.2
    pointsize = radius / 8.0

    # first disc
    println("--------- First disc ---------")
    center = [0.2, 0.2]
    firstdisc = creatediscmaterialdomain(center, radius, pointsize, [0.1, 0.1], g, ρ, E, ν)

    # second disc
    println("--------- Second disc ---------")
    center = [0.8, 0.8]
    seconddisc =
        creatediscmaterialdomain(center, radius, pointsize, [-0.1, -0.1], g, ρ, E, ν)

    # all material points
    materialpoints = [firstdisc seconddisc]
    totalmass = 0.0
    for point in materialpoints
        totalmass += point.m
    end
    println("--------- Total material points ---------")
    println("  mass: ", totalmass)
    println("  number of material points: ", length(materialpoints))
    println()

    # time stepping loop
    dt = 1.0e-3
    t_f = 3.5e0
    println("--------- Time stepping loop ---------")
    println("  start time: ", 0.0)
    println("  stop time: ", t_f)
    println("  step size: ", dt)
    println()
    for t = 0.0:dt:t_f
        # reset grid
        resetgrid(grid)

        # material points to grid
        for materialpoint in materialpoints
            transfermaterialpointtogrid(materialpoint, grid)
        end

        # solve grid momentum
        for gridpoint in grid.points
            gridpoint.p += gridpoint.f * dt

            # boundary conditions
            if (gridpoint.isfixed[1])
                gridpoint.p[1] = 0.0
                gridpoint.f[1] = 0.0
            end
            if (gridpoint.isfixed[2])
                gridpoint.p[2] = 0.0
                gridpoint.f[2] = 0.0
            end
        end

        # grid to material points
        for materialpoint in materialpoints
            transfergridtomaterialpoint(materialpoint, dt, grid)
        end
    end
end

# ------------------------------------------------------------------------------
# run main
# ------------------------------------------------------------------------------

main()

# ------------------------------------------------------------------------------
