# ------------------------------------------------------------------------------
# Bases
# ------------------------------------------------------------------------------

function getbasismatrices(materialpoint::MaterialPoint, gridpoint::GridPoint, grid::Grid)
    distance = materialpoint.x - gridpoint.x

    # interpolation
    shape = zeros(2)
    shape[1] = 1.0 - abs(distance[1]) / grid.dx[1]
    shape[2] = 1.0 - abs(distance[2]) / grid.dx[2]

    # gradient
    gradient = zeros(2)
    gradient[1] = -shape[2] * sign(distance[1]) / grid.dx[1]
    gradient[2] = -shape[1] * sign(distance[2]) / grid.dx[2]

    # return
    (shape[1] * shape[2], gradient)
end

# ------------------------------------------------------------------------------
