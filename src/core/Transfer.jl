# ------------------------------------------------------------------------------
# Transfer between material points and background grid
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# utility
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

# ------------------------------------------------------------------------------
# material point to grid transfer
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# grid to material point transfer
# ------------------------------------------------------------------------------

function transfergridtomaterialpoint(materialpoint::MaterialPoint, dt::Float64, grid::Grid)
    # transfer to material point from all adjacent grid points
    adjacentgridindices = getadjacentgridindices(materialpoint, grid)

    # accumulate adjustments from grid points
    dx = zeros(2)
    for index in adjacentgridindices
        gridpoint = grid.points[index]

        interpolation, gradient = getbasismatrices(materialpoint, gridpoint, grid)

        v = zeros(2)
        if (gridpoint.m > eps(Float64) * 100.0)
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
# grid to material point transfer with yield
# ------------------------------------------------------------------------------

function transfergridtomaterialpointwithyieldpass1(
    materialpoint::MaterialPoint,
    dt::Float64,
    grid::Grid,
)
    # transfer to material point from all adjacent grid points
    adjacentgridindices = getadjacentgridindices(materialpoint, grid)

    # accumulate adjustments from grid points
    for index in adjacentgridindices
        gridpoint = grid.points[index]

        if (gridpoint.m > eps(Float64) * 100.0)
            interpolation, gradient = getbasismatrices(materialpoint, gridpoint, grid)
            materialpoint.v += dt * (interpolation * gridpoint.f / gridpoint.m)
        end
    end

    # adjust grid point velocity
    for index in adjacentgridindices
        gridpoint = grid.points[index]

        interpolation, gradient = getbasismatrices(materialpoint, gridpoint, grid)
        gridpoint.v += interpolation * materialpoint.m * materialpoint.v / gridpoint.m
    end
end

function transfergridtomaterialpointwithyieldpass2(
    materialpoint::MaterialPoint,
    dt::Float64,
    grid::Grid,
)
    # transfer to material point from all adjacent grid points
    adjacentgridindices = getadjacentgridindices(materialpoint, grid)

    # accumulate adjustments from grid points
    dx = zeros(2)
    for index in adjacentgridindices
        gridpoint = grid.points[index]

        interpolation, gradient = getbasismatrices(materialpoint, gridpoint, grid)

        if (gridpoint.m > eps(Float64) * 100.0)
            dx += dt * (interpolation * gridpoint.p / gridpoint.m)
            materialpoint.dF += dt * gridpoint.v * gradient'
        end
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
    dσ, dε_plastic, dα = incrementstrainplastic(materialpoint)
    ε_current = materialpoint.ε
    σ_current = materialpoint.σ
    ε_plastic_current = materialpoint.ε_plastic
    materialpoint.σ += dσ
    materialpoint.ε_plastic += dε_plastic
    materialpoint.α += dα

    # volume and momentum
    materialpoint.V = det(materialpoint.F) * materialpoint.V_0
    materialpoint.p = materialpoint.v * materialpoint.m
end

# ------------------------------------------------------------------------------
