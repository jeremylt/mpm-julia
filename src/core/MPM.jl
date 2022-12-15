# ------------------------------------------------------------------------------
# Collects all files needed for main drivers
# ------------------------------------------------------------------------------

# logging level
quiet = false

# material points
include("MaterialPoint.jl")
include("BoxMaterialDomain.jl")
include("DiscMaterialDomain.jl")

# background grid
include("GridPoint.jl")
include("Grid.jl")

# material point and grid transfer
include("Basis.jl")
include("Transfer.jl")

# ------------------------------------------------------------------------------
