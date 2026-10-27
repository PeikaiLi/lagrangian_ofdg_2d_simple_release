module DGOperator

import ..Basis
import ..Euler
import ..FaceFlux
import ..Geometry
import ..Problems
import ..Discretization
import ..MeshMotion
import ..Meshes
import ..States

# All arrays describe the same RK stage; time is the stage time used by the source.
# Return a new basis_count-by-4-by-cell_count volume contribution to dQ/dt.
function assemble_volume_rhs(
    state::States.DGState,
    case::Discretization.DGCase,
    time::Real,
    cells_geometries::Vector{Geometry.CellGeometry},
    cells_coefficients::Array{Float64,3},
    vertices_velocities_mat::Matrix{Float64},
)

    mesh = case.mesh
    problem = case.problem
    options = case.options
    cells_gamma_ls = case.cells_gamma_ls

    cells_volume_residual_array = zeros(size(state.cells_moments_array))
    exponents_ls = Basis.basis_powers(options.degree)
    
end