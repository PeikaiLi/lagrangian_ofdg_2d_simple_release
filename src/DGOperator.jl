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

    # No source for this version
    function assemble_volume_rhs(
        state::States.DGState,
        case::Discretization.DGCase,
        time::Real, # for future source version
        cells_geometries::Vector{Geometry.CellGeometry},
        cells_coefficients::Array{Float64,3},
        vertices_velocities_mat::Matrix{Float64},
    )

        num_of_basis, number_of_cons, num_of_cell = size(state.cells_moments_array)
        cells_volume_rhs_array = zeros(num_of_basis, number_of_cons, num_of_cell)
        cells_gamma_ls = case.cells_gamma_ls
        exponents_ls = Basis.basis_powers(case.options.degree)

        for k in 1:num_of_cell

            gamma = cells_gamma_ls[k]
            geometry = cells_geometries[k]
            coefficients_mat = cells_coefficients[:,:,k]
            cell_velocities_mat = vertices_velocities_mat[:,case.mesh.cells_vertex_indices_mat[:, k]]
            num_of_quadrature_pts = length(geometry.quadrature_weights_ls)

        
            cell_motion_cache =  MeshMotion.build_CellMotionCache(     
                geometry,
                cell_velocities_mat,
                case.options,
            )

            for pt_index in 1:num_of_quadrature_pts
                pt_coord = geometry.quadrature_coords_ls[pt_index]
                cons_vec = Geometry.eval_cons(coefficients_mat, geometry, pt_coord, exponents_ls)
                rho, velocity_x, velocity_y, pressure_value = Euler.cons2prim(cons_vec, gamma)
                # Pressure fluxes Gx and Gy in the paper.
                pressure_flux_x_vec = [
                    0.0,
                    pressure_value,
                    0.0,
                    pressure_value * velocity_x,
                ]

                pressure_flux_y_vec = [
                    0.0,
                    0.0,
                    pressure_value,
                    pressure_value * velocity_y,
                ]

                temp = 
                    velocity_x*geometry.basis_dx_at_quadrature_mat[:,pt_index] + 
                    velocity_y*geometry.basis_dy_at_quadrature_mat[:,pt_index] + 
                    cell_motion_cache.dt_basis_at_quadrature_mat[:,pt_index]

                for i in 1:number_of_cons
                    cells_volume_rhs_array[:,i,k] += (
                        cons_vec[i] * temp + 
                        pressure_flux_x_vec[i] * geometry.basis_dx_at_quadrature_mat[:,pt_index] + 
                        pressure_flux_y_vec[i] * geometry.basis_dy_at_quadrature_mat[:,pt_index]
                    ) * geometry.quadrature_weights_ls[pt_index]
                end

            end

        end 
        
        return cells_volume_rhs_array
    end
    
end # module DGOperator