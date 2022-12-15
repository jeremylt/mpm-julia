# ------------------------------------------------------------------------------
# 2D example - fixed beam
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
#   Note: Plots needs to be installed to use
# ------------------------------------------------------------------------------

using Printf
using Plots

# ------------------------------------------------------------------------------
# object definitions
# ------------------------------------------------------------------------------

include("core/MPM.jl")

# ------------------------------------------------------------------------------
# driver
# ------------------------------------------------------------------------------

function main()
    # physical constants
    g = -1.0
    ρ = 1000.0 # density
    E = 1000.0 # Young's modulus
    ν = 0.3    # Poisson's ratio

    # setup grid
    println("--------- Background grid ---------")
    numcells = 20
    grid = Grid(1.0, 1.0, numcells + 1, numcells + 1)

    # boundary
    for gridpoint in grid.points
        if (gridpoint.x[1] == 0.0)
            gridpoint.isfixed[1] = true
            gridpoint.isfixed[2] = true
        end
    end

    # setup disc parameters
    width = 0.5
    height = 0.25
    pointsize = height / 8.0

    # beam
    println("--------- Beam ---------")
    lowerleftcorner = [0.0, 0.5]
    materialpoints = createboxmaterialdomain(
        lowerleftcorner,
        width,
        height,
        pointsize,
        [0.0, 0.0],
        g,
        ρ,
        E,
        ν,
    )

    # all material points
    totalmass = 0.0
    for point in materialpoints
        totalmass += point.m
    end
    println("--------- Total material points ---------")
    println("  mass: ", totalmass)
    println("  number of material points: ", length(materialpoints))
    println()

    # plotted quantities
    plotincrement = 1
    times = []
    strainenergies = []
    kineticenergies = []

    # time stepping loop
    dt = 1.0e-3
    t_f = 3.8e0
    println("--------- Time stepping loop ---------")
    println("  start time: ", 0.0)
    println("  stop time: ", t_f)
    println("  step size: ", dt)
    println()
    step = 0
    for t = 0.0:dt:t_f
        println("  current time: ", t)

        # plot current locations
        if (step % plotincrement == 0)
            println("  plotting current location")
            materialpoints_x = Array{Float64}(undef, length(materialpoints))
            materialpoints_y = Array{Float64}(undef, length(materialpoints))
            for i = 1:length(materialpoints)
                materialpoints_x[i] = materialpoints[i].x[1]
                materialpoints_y[i] = materialpoints[i].x[2]
            end
            scatter(
                materialpoints_x,
                materialpoints_y,
                aspect_ratio = :equal,
                lims = [0.0, 1.0],
                title = "Two Disc MPM",
                label = "",
                xlabel = "x",
                ylabel = "y",
            )
            savefig("Beam-$(@sprintf("%0.3f", t)).png")
        end

        @time "  material points to grid points" begin
            # reset grid
            resetgrid(grid)

            # material points to grid
            for materialpoint in materialpoints
                transfermaterialpointtogrid(materialpoint, grid)
            end
        end

        # solve grid momentum
        @time "  grid point solve" begin
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
        end

        # grid to material points
        @time "  grid points to material points" begin
            for materialpoint in materialpoints
                transfergridtomaterialpoint(materialpoint, dt, grid)
            end
        end

        # update energy plot data

        if (step % plotincrement == 0)
            strainenergy = 0.0
            kineticenergy = 0.0
            for materialpoint in materialpoints
                for i = 1:3
                    strainenergy +=
                        0.5 * materialpoint.σ[i] * materialpoint.ε[i] * materialpoint.V
                end
                kineticenergy +=
                    0.5 * (materialpoint.v[1]^2 + materialpoint.v[2]^2) * materialpoint.m
            end
            push!(times, t)
            push!(strainenergies, strainenergy)
            push!(kineticenergies, kineticenergy)
        end

        # log progress
        p = zeros(2)
        for materialpoint in materialpoints
            p += materialpoint.p
        end
        totalmass = 0.0
        for point in materialpoints
            totalmass += point.m
        end
        println("  momentum (x, y): (", p[1], ", ", p[2], ")")
        println("  mass: ", totalmass)
        println()

        # increment plot step counter
        step += 1
    end

    # plotting
    plot(
        times,
        [strainenergies, kineticenergies],
        title = "Strain Energy and Kinetic Energy",
        label = ["Strain Energy" "Kinetic Energy"],
        xlabel = "time",
    )
    savefig("Beam-StrainEnergyAndKineticEnergy.png")
end

# ------------------------------------------------------------------------------
# run main
# ------------------------------------------------------------------------------

main()

# ------------------------------------------------------------------------------
