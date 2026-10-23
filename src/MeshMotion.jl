module MeshMotion

    using LinearAlgebra: dot, norm
    import ..Basis
    import ..Euler
    import ..FaceFlux
    import ..Meshes
    import ..Discretization
    import ..Geometry
    import ..States

    # Use the current stage state geometry coefficient...
    # Return velocity matrix; 2*number_of_vertices
    function compute_vertex_velocities(
        state::States.DGState,
        case::Discretization.DGCase,
        cells_geometries::Vector{Geometry.CellGeometry},
        cells_coefficients::Array{Float64,3},
    )

        mesh = case.mesh
        options = case.options

        exponents_ls = Basis.basis_powers(options.degree)
        number_of_velocity_groups = maximum(mesh.vertex2group_index)

        group_velocity_sums_mat = zeros(2, number_of_velocity_groups)
        group_weight_sums_ls = zeros(number_of_velocity_groups)
        wall_group_flags_ls = falses(number_of_velocity_groups)

        for face in mesh.faces_ls
            if face.right_cell_index == 0 && face.boundary == :wall
                start_group_index = mesh.vertex2group_index[face.start_vertex_index]
                end_group_index = mesh.vertex2group_index[face.end_vertex_index]
                wall_group_flags_ls[start_group_index] = true
                wall_group_flags_ls[end_group_index] = true
            end
        end

        for face in mesh.faces_ls
            _, _, face_length, normal = Meshes.current_face_geometry(state.vertices_coords_mat, face)
            if face_length == 0.0
                continue
            end

            tangent = [-normal[2], normal[1]]

            for vertex_index in (face.start_vertex_index, face.end_vertex_index)
                vertex_coord = state.vertices_coords_mat[:, vertex_index]

                left_cons_vec,
                right_cons_vec,
                left_gamma,
                right_gamma = FaceFlux.face_cons(
                    face,
                    vertex_coord,
                    cells_coefficients,
                    cells_geometries,
                    case.problem,
                    case.cells_gamma_ls,
                    exponents_ls,
                    normal,
                )

                star_normal_velocity, _ = Euler.hllc_star(left_cons_vec, right_cons_vec, left_gamma, right_gamma, normal)

                left_velocity = left_cons_vec[2:3] / left_cons_vec[1]
                right_velocity = right_cons_vec[2:3] / right_cons_vec[1]

                tangent_velocity = dot(0.5 * (left_velocity + right_velocity), tangent)

                candidate_velocity = star_normal_velocity * normal + tangent_velocity * tangent

                group_index = mesh.vertex2group_index[vertex_index]
                weight = 1.0
                if wall_group_flags_ls[group_index] && face.right_cell_index > 0
                    weight = 2.0
                end

                group_velocity_sums_mat[:, group_index] += weight * candidate_velocity
                group_weight_sums_ls[group_index] += weight
            end
        end

        vertices_velocities_mat = zeros(size(state.vertices_coords_mat))
        for vertex_index in 1:size(state.vertices_coords_mat, 2)
            group_index = mesh.vertex2group_index[vertex_index]
            if group_weight_sums_ls[group_index] > 0.0
                vertices_velocities_mat[:, vertex_index] = group_velocity_sums_mat[:, group_index] /
                    group_weight_sums_ls[group_index]
            end
        end

        # Remove the normal mesh velocity at stationary walls; keep tangential motion.
        for face in mesh.faces_ls
            if face.right_cell_index == 0 && face.boundary == :wall
                _, _, face_length, normal = Meshes.current_face_geometry(state.vertices_coords_mat, face)
                if face_length == 0.0
                    continue
                end

                for vertex_index in (face.start_vertex_index, face.end_vertex_index)
                    vertices_velocities_mat[:, vertex_index] -=
                        dot(vertices_velocities_mat[:, vertex_index], normal) * normal
                end
            end
        end

        return vertices_velocities_mat # for all vex
    end

    # basis_time_derivative at a point in a given cell
    function basis_time_derivative_at_one_pt(
        geometry::Geometry.CellGeometry,
        cell_velocities_mat::Matrix{Float64}, # counterclock
        exponents_ls::Vector{Tuple{Int,Int}},
        pt_coord::Vector{Float64},
        quadrature_order
    )
        quadrature1d_pts_ls, quadrature1d_wts_ls = gauss_quad0_1(quadrature_order)
        # boundary integral for each edge
        # I should move boundary integral quadraturepts to Geometry first other write it will repeat too much times
        basis_dt_vec = zeros(length(exponents_ls))

        raw_basis_val_at_pt = eval_basis_values(geometry, pt_coord, exponents_ls) + geometry.raw_basis_means_ls

        dt_area = 0.0
        boundary_integrals = zeros(length(exponents_ls))
        



        return basis_dt_mat
    end

end # module MeshMotion
