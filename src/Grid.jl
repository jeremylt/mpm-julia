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
    v_cell::Array{Float64} # size of each cell/element
    v_cell_inv::Array{Float64} # inverse of size of each cell
    points::Array{GridPoint} # array of all grid points

    # constructor
    Grid(length_x::Float64, length_y::Float64, numnodes_x::Int64, numnodes_y::Int64) = (
        # cell sizes
        v_cell = zeros(2);
        v_cell_inv = zeros(2);
        v_cell[1] = length_x / Float64(numnodes_x - 1.0);
        v_cell[2] = length_y / Float64(numnodes_y - 1.0);
        v_cell_inv[1] = 1.0 / v_cell[1];
        v_cell_inv[2] = 1.0 / v_cell[2];

        # grid points array
        numnodes = numnodes_x * numnodes_y;
        gridpoints = Array{GridPoint}(undef, numnodes);
        for j = 1:numnodes_y
            y = (j - 1) * v_cell[2]
            for i = 1:numnodes_x
                x = (i - 1) * v_cell[1]
                index = index2Dto1D(i, j, numnodes_x, numnodes_y)
                gridpoints[index] = GridPoint(x, y)
            end
        end;

        # constructor
        new(
            [length_x, length_y],
            [numnodes_x, numnodes_y],
            numnodes,
            v_cell,
            v_cell_inv,
            gridpoints,
        )
    )
end

# ------------------------------------------------------------------------------
# material point to grid transfer
# ------------------------------------------------------------------------------

function getadjacentgridindices(materialpoint::MaterialPoint, grid::Grid)
    # cell sizes
    celllength_x = grid.v_cell[1]
    celllength_y = grid.v_cell[2]

    # find row/column index
    cellcolumn = floor(materialpoint.x[1] * celllength_x) + 1
    cellrow = floor(materialpoint.x[2] * celllength_y) + 1

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

# ------------------------------------------------------------------------------
