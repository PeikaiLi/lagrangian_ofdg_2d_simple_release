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
function assemble_volume_residual(
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

    for cell_index in eachindex(cells_geometries)
        geometry = cells_geometries[cell_index]
        cell_vertex_indices_ls = mesh.cells_vertex_indices_mat[:, cell_index]
        cell_velocities_mat = vertices_velocities_mat[:, cell_vertex_indices_ls]

        # Geometry and motion caches belong to this cell at the current RK stage.
        motioncache = MeshMotion.build_CellMotionCache(
            geometry, cell_velocities_mat, options,
        )
        basis_dt_mat = motioncache.dt_basis_at_quadrature_mat

        for quadrature_point_index in eachindex(geometry.quadrature_weights_ls)
            cons_vec = transpose(cells_coefficients[:, :, cell_index]) *
                geometry.basis_at_quadrature_mat[:, quadrature_point_index]

            flux_x_ls, flux_y_ls = Euler.euler_flux(
                cons_vec,
                cells_gamma_ls[cell_index],
            )

            source_vec = Problems.source_term(
                problem,
                geometry.quadrature_coords_ls[quadrature_point_index],
                time,
            )

            for basis_index in eachindex(exponents_ls)
                basis_value = geometry.basis_at_quadrature_mat[basis_index, quadrature_point_index]
                basis_dx = geometry.basis_dx_at_quadrature_mat[basis_index, quadrature_point_index]
                basis_dy = geometry.basis_dy_at_quadrature_mat[basis_index, quadrature_point_index]
                basis_dt = basis_dt_mat[basis_index, quadrature_point_index]
                integrand = basis_dx * flux_x_ls + basis_dy * flux_y_ls +
                    basis_dt * cons_vec + basis_value * source_vec

                cells_volume_residual_array[basis_index, :, cell_index] +=
                    geometry.quadrature_weights_ls[quadrature_point_index] * integrand
            end
        end
    end

    return cells_volume_residual_array
end
