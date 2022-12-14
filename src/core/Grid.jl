# ------------------------------------------------------------------------------
# Background grid
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# utility
# ------------------------------------------------------------------------------

function index2Dto1D(i::Int64, j::Int64, numcolumns::Int64, numrows::Int64)
    # convert to 1D index
    index = numcolumns * (j - 1) + i

    # check bounds
    if (index < 1 || index > numrows * numcolumns)
        throw(
            DomainError(
                i,
                j,
                index,
                numrows * numcolumns,
                "Index does not fall within 1D range",
            ),
        )
    end

    # return
    index
end

# ------------------------------------------------------------------------------
# grid
# ------------------------------------------------------------------------------

mutable struct Grid
    # internal state data
    length::Array{Float64} # length in x and y directions
    shape::Array{Int64}
    numnodes::Int64
    dx::Array{Float64} # size of each cell/element
    dx_inv::Array{Float64} # inverse of size of each cell
    points::Array{GridPoint} # array of all grid points

    # constructor
    Grid(length_x::Float64, length_y::Float64, numnodes_x::Int64, numnodes_y::Int64) = (
        # cell sizes
        dx = zeros(2);
        dx_inv = zeros(2);
        dx[1] = length_x / Float64(numnodes_x - 1.0);
        dx[2] = length_y / Float64(numnodes_y - 1.0);
        dx_inv[1] = 1.0 / dx[1];
        dx_inv[2] = 1.0 / dx[2];

        # grid points array
        numnodes = numnodes_x * numnodes_y;
        gridpoints = Array{GridPoint}(undef, numnodes);
        for j = 1:numnodes_y
            y = (j - 1) * dx[2]
            for i = 1:numnodes_x
                x = (i - 1) * dx[1]
                index = index2Dto1D(i, j, numnodes_x, numnodes_y)
                gridpoints[index] = GridPoint(x, y)
            end
        end;

        # print information
        println("Background grid:");
        println("  x range: ", 0.0, " - ", length_x);
        println("  y range: ", 0.0, " - ", length_y);
        println("  number columns (x): ", numnodes_x);
        println("  number rows (y): ", numnodes_y);
        println("  total gridpoints: ", numnodes);
        println();

        # constructor
        new(
            [length_x, length_y],
            [numnodes_x, numnodes_y],
            numnodes,
            dx,
            dx_inv,
            gridpoints,
        )
    )
end

# ------------------------------------------------------------------------------
# reset grid
# ------------------------------------------------------------------------------

function resetgrid(grid::Grid)
    for point in grid.points
        point.m = 0.0
        point.p = [0.0, 0.0]
        point.f = [0.0, 0.0]
    end
end

# ------------------------------------------------------------------------------
# material point to grid transfer
# ------------------------------------------------------------------------------

function getadjacentgridindices(materialpoint::MaterialPoint, grid::Grid)
    # find row/column index
    cellcolumn = Int64(floor(materialpoint.x[1] * grid.dx_inv[1]) + 1)
    cellrow = Int64(floor(materialpoint.x[2] * grid.dx_inv[2]) + 1)

    # check bounds
    if (cellcolumn < 1 || cellcolumn > grid.shape[1])
        throw(DomainError(cellcolumn, "Material point outside of the grid"))
    end
    if (cellrow < 1 || cellrow > grid.shape[2])
        throw(DomainError(cellrow, "Material point outside of the grid"))
    end

    # convert to 1D index
    firstpoint = index2Dto1D(cellcolumn, cellrow, grid.shape[1], grid.shape[2])
    secondpoint = index2Dto1D(cellcolumn, cellrow + 1, grid.shape[1], grid.shape[2])

    # return
    [firstpoint, secondpoint, firstpoint + 1, secondpoint + 1]
end

function transfermaterialpointtogrid(materialpoint::MaterialPoint, grid::Grid)
    # transfer material point to all adjacent grid points
    adjacentgridindices = getadjacentgridindices(materialpoint, grid)
    for index in adjacentgridindices
        gridpoint = grid.points[index]

        interpolation, gradient = getbasismatrices(materialpoint, gridpoint, grid)

        gridpoint.m += interpolation * materialpoint.m
        gridpoint.p += interpolation * materialpoint.m * materialpoint.v
        gridpoint.f += materialpoint.f_external
        gridpoint.f[1] -=
            materialpoint.V *
            (gradient[1] * materialpoint.σ[1] + gradient[2] * materialpoint.σ[3])
        gridpoint.f[2] -=
            materialpoint.V *
            (gradient[2] * materialpoint.σ[2] + gradient[1] * materialpoint.σ[3])
    end
end

function transfergridtomaterialpoint(materialpoint::MaterialPoint, dt::Float64, grid::Grid)
    # transfer to material point from all adjacent grid points
    adjacentgridindices = getadjacentgridindices(materialpoint, grid)

    # accumulate adjustments from grid points
    dx = zeros(2)
    for index in adjacentgridindices
        gridpoint = grid.points[index]

        interpolation, gradient = getbasismatrices(materialpoint, gridpoint, grid)

        v = zeros(2)
        if (gridpoint.m > 1.0e-8)
            v = gridpoint.p / gridpoint.m
            materialpoint.v += dt * (interpolation * gridpoint.f / gridpoint.m)
            dx += dt * (interpolation * gridpoint.p / gridpoint.m)
        end
        materialpoint.dF += dt * v * gradient'
    end

    # update material point
    materialpoint.x += dx
    materialpoint.F = materialpoint.dF * materialpoint.F

    # stress
    dε = zeros(3)
    dε[1] = materialpoint.dF[1, 1] - 1.0
    dε[2] = materialpoint.dF[2, 2] - 1.0
    dε[3] = materialpoint.dF[1, 2] + materialpoint.dF[2, 1]
    materialpoint.ε += dε
    materialpoint.dF = I(2)

    # strain
    E = materialpoint.E
    ν = materialpoint.ν
    k = E / (1.0 + ν) / (1.0 - 2.0 * ν)
    materialpoint.σ[1] += k * ((1.0 - ν) * dε[1] + ν * dε[2])
    materialpoint.σ[2] += k * ((1.0 - ν) * dε[2] + ν * dε[1])
    materialpoint.σ[3] += k * ((0.5 - ν) * dε[3])

    # volume and momentum
    materialpoint.V = det(materialpoint.F) * materialpoint.V_0
    materialpoint.p = materialpoint.v * materialpoint.m
end

# ------------------------------------------------------------------------------
