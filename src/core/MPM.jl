# ------------------------------------------------------------------------------
# Collects all files needed for main drivers
# ------------------------------------------------------------------------------

# Material points
include("MaterialPoint.jl")
include("BoxMaterialDomain.jl")
include("DiscMaterialDomain.jl")

# Background grid
include("GridPoint.jl")
include("Grid.jl")

# Material point and grid transfer
include("Basis.jl")
include("Transfer.jl")

# ------------------------------------------------------------------------------
