# Load all 
include("Meshes.jl")
include("Quadrature.jl")
include("Options.jl")
include("Basis.jl")
include("Geometry.jl")
include("States.jl")
include("Euler.jl")
include("Problems.jl")
include("Discretization.jl")
include("FaceFlux.jl")
include("MeshMotion.jl")
include("DGOperator.jl")
include("Positivity.jl")

include("OF.jl")
include("TimeIntegration.jl")
include("Solver.jl")
include("Diagnostics.jl")
include("Visualization.jl") 

include("cases/ConstantCase.jl")
include("cases/IsentropicCase.jl")
include("cases/LaxCase.jl") 
include("cases/SodCase.jl") 
include("cases/QuadrantCase.jl") 


include("cases/BlastWaveCase.jl")
include("cases/ShuOsherCase.jl")

println("Loaded successfully.")

