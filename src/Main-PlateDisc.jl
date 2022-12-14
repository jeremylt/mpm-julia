# ------------------------------------------------------------------------------
# 2D example - plate and disc
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
#   Note: needs to be installed to use
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
    g = 0.0

    # setup grid
    println("--------- Background grid ---------")
    numcells = 25
    grid = Grid(60.0, 60.0, numcells + 1, numcells + 1)

    # boundary
    for gridpoint in grid.points
        if (gridpoint.x[1] == 0.0 || gridpoint.x[1] == 60.0 || gridpoint.x[2] == 0.0)
            gridpoint.isfixed[1] = true
            gridpoint.isfixed[2] = true
        end
    end

    # Plate
    println("--------- Plate ---------")
    ρ = 2700e-12    # density
    E = 78.2e3      # Young's modulus
    ν = 0.3         # Poisson's ratio
    σ_yield = 300.0 # yield stress
    lowerleftcorner = [0.0, 0.0]
    width = 60.0
    height = 40.0
    pointsize = width / 50
    plate = createboxmaterialdomain(
        lowerleftcorner,
        width,
        height,
        pointsize,
        [0.0, 0.0],
        g,
        ρ,
        E,
        ν,
        σ_yield,
    )

    # Disc
    ρ = 7850e-12     # density
    E = 200.0e3      # Young's modulus
    ν = 0.3          # Poisson's ratio
    σ_yield = 1.0e24 # yield stress
    radius = 9.6 / 2.0
    println("--------- Disc ---------")
    center = [30.0, 50.0]
    disc = creatediscmaterialdomain(
        center,
        radius,
        pointsize,
        [0.0, -1160.0e3],
        g,
        ρ,
        E,
        ν,
        σ_yield,
    )

    # all material points
    materialpoints = [plate; disc]
    totalmass = 0.0
    for point in materialpoints
        totalmass += point.m
    end
    println("--------- Total material points ---------")
    println("  mass: ", totalmass)
    println("  number of material points: ", length(materialpoints))
    println()

    # plotted quantities
    plotincrement = 100
    times = []
    strainenergies = []
    kineticenergies = []

    # time stepping loop
    dt = 1.0e-8
    t_f = 1.0e-4
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
            materialpoints_x = [materialpoint.x[1] for materialpoint in materialpoints]
            materialpoints_y = [materialpoint.x[2] for materialpoint in materialpoints]
            scatter(
                materialpoints_x,
                materialpoints_y,
                aspect_ratio = :equal,
                lims = [0.0, 60.0],
                title = "Plate and Disc MPM",
                label = "",
                xlabel = "x",
                ylabel = "y",
            )
            savefig("PlateDisc-$(@sprintf("%0.8f", t)).png")
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
            # first pass
            for materialpoint in materialpoints
                transfergridtomaterialpointwithyieldpass1(materialpoint, dt, grid)
            end

            # boundary conditions
            for gridpoint in grid.points
                if (gridpoint.isfixed[1])
                    gridpoint.v[1] = 0.0
                end
                if (gridpoint.isfixed[2])
                    gridpoint.v[2] = 0.0
                end
            end

            # second pass
            for materialpoint in materialpoints
                transfergridtomaterialpointwithyieldpass2(materialpoint, dt, grid)
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
    savefig("PlateDisc-StrainEnergyAndKineticEnergy.png")
end

# ------------------------------------------------------------------------------
# run main
# ------------------------------------------------------------------------------

main()

# ------------------------------------------------------------------------------