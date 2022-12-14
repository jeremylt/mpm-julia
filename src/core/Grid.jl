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
        point.v = [0.0, 0.0]
        point.p = [0.0, 0.0]
        point.f = [0.0, 0.0]
    end
end

# ------------------------------------------------------------------------------
