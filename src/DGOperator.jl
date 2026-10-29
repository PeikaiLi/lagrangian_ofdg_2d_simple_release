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

        
            cell_motion_cache = MeshMotion.build_CellMotionCache(     
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


    function eval_a_face_flux(
        face::Meshes.Face,
        state::States.DGState,
        case::Discretization.DGCase,
        problem::Problems.EulerProblem,
        cells_geometries::Vector{Geometry.CellGeometry},
        cells_coefficients::Array{Float64,3},
        vertices_velocities_mat::Matrix{Float64},
        )
       

        start_coord, 
        end_coord, 
        edge_length, 
        normal = 
            current_face_geometry(state.vertices_coords_mat, face)

        cells_gamma_ls = case.cells_gamma_ls

        exponents_ls = Basis.basis_powers(case.options.degree)

        edge_quadrature_weights_ls = case.options.degree.quadrature1d_wts_ls
        num_of_quadrature_pt = length(edge_quadrature_weights_ls)

        edge_quadrature_weights_ls = edge_quadrature_weights_ls *edge_length

        edge_quadrature_coords_mat = zeros(2,num_of_quadrature_pt)
        advection_flux_mat = zeros(4,num_of_quadrature_pt)
        pressure_flux_mat = zeros(4,num_of_quadrature_pt)
        integrated_pressure_jump = zeros(num_of_quadrature_pt)
        upwind_pressure_max = 0

        for pt_index in 1:num_of_quadrature_pt
            point_coord = case.options.degree.quadrature1d_pts_ls[:,pt_index]
            edge_quadrature_coords_mat[:,pt_index] = point_coord
            left_cons_vec, right_cons_vec, left_gamma, right_gamma = 
                face_cons(face, point_coord, cells_coefficients, cells_geometries,
                    problem, cells_gamma_ls, exponents_ls, normal)
            
            advection_flux, pressure_flux, oriented_pressure_jump, upwind_pressure = 
                moving_face_flux(
                    left_cons_vec, 
                    right_cons_vec, 
                    left_gamma, 
                    right_gamma,
                    normal, 
                    mesh_velocity
                )

            advection_flux_mat[:, pt_index] = advection_flux
            pressure_flux_mat[:, pt_index] = pressure_flux
            integrated_pressure_jump += edge_quadrature_weights_ls[pt_index] * oriented_pressure_jump
            upwind_pressure_max = max(upwind_pressure_max, upwind_pressure)
        end

        #KXRCF 
        is_shock_face = false
        if !case.problem.smooth && case.options.degree >= 2
            detector_scale = edge_length^((case.options.degree + 3) / 2.0) *
                max(upwind_pressure_max, 1.0e-14)
            is_shock_face = abs(integrated_pressure_jump) > detector_scale
        end

        return (
            edge_quadrature_coords_mat = quadrature_coords_mat,
            edge_quadrature_weights_ls = edge_quadrature_weights_ls,
            advection_flux_mat = advection_flux_mat,
            pressure_flux_mat = pressure_flux_mat,
            is_shock_face = is_shock_face,
        )
    end

    

    function eval_all_face_flux(
        # cells_volume_rhs_array
        state::States.DGState,
        case::Discretization.DGCase,
        problem::Problems.EulerProblem,
        cells_geometries::Vector{Geometry.CellGeometry},
        cells_coefficients::Array{Float64,3},
        vertices_velocities_mat::Matrix{Float64},
    )

        num_of_face = length(case.mesh.faces_ls)
        num_of_quadrature_pt = length(case.options.degree.quadrature1d_wts_ls)


        for face_index in 1:num_of_face

            face_data = eval_a_face_flux(
                case.mesh.faces_ls[face_index],
                state,
                case,
                problem,
                cells_geometries,
                cells_coefficients,
                vertices_velocities_mat,
            )



            for pt_index in 1:num_of_quadrature_pt


            
            
            end
        end

        
    end
    


end # module DGOperator