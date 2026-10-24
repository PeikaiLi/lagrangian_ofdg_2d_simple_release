module MeshMotion

    using LinearAlgebra: dot, norm
    import ..Options
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
        cell_velocities_mat::Matrix{Float64}, # same order and size as cell_vertices_coords_mat
        exponents_ls::Vector{Tuple{Int,Int}},
        pt_coord::Vector{Float64},
        options::Options.DGOptions
    )
               
        basis_dt_vec_at_pt = zeros(length(exponents_ls))
        basis_values_at_pt = eval_basis_values(geometry, pt_coord, exponents_ls) 
        raw_basis_val_at_pt = basis_values_at_pt + geometry.raw_basis_means_ls
        xi = raw_basis_val_at_pt[2]
        eta = raw_basis_val_at_pt[3]

        
        
        # which is a boundary integral ug*n
        num_of_edge = size(geometry.cell_vertices_coords_mat)[2]
        # I use this to deal with future triangle mesh

        quadrature1d_pts_ls,quadrature1d_wts_ls =  options.quadrature1d_pts_ls, options.quadrature1d_wts_ls

        quadrature_coords_at_edge = zeros(2,num_of_edge*length(quadrature1d_pts_ls))
        quadrature_weight_at_edge = zeros(num_of_edge*length(quadrature1d_pts_ls))
        u_g_quadrature_at_edge = zeros(2,num_of_edge*length(quadrature1d_pts_ls))
        norm_quadrature_at_edge = zeros(2,num_of_edge*length(quadrature1d_pts_ls))

        pt_index = 0
        for start_vex_index in 1:num_of_edge
            pt_index += 1
            start_vex_coord = geometry.cell_vertices_coords_mat[:,start_vex_index]
            end_vex_index = mod(start_vex_index+1,num_of_edge)
            end_vex_coord = geometry.cell_vertices_coords_mat[:,end_vex_index]

            start_vex_volcity = cell_velocities_mat[:,start_vex_index]
            end_vex_volcity = cell_velocities_mat[:,end_vex_index]

            tangent_vector = end_vex_coord - start_vex_coord
            length_tangent_vector = norm(tangent_vector)
            Options.check_positive(length_tangent_vector, "zero edge in moving basis")
            tangent_vector = tangent_vector/norm(tangent_vector) 
            normal_vector = [-tangent_vector[2],tangent_vector[1]] # clock 90 degree

            for i in 1:length(quadrature1d_pts_ls)
                quadrature_coords_at_edge[:,pt_index] = start_vex_coord + quadrature1d_pts_ls[i](end_vex_coord - start_vex_coord)
                quadrature_weight_at_edge[:,pt_index] = quadrature1d_wts_ls[i]
                u_g_quadrature_at_edge[:,pt_index] = start_vex_volcity + quadrature1d_pts_ls[i](end_vex_volcity - start_vex_volcity)
                norm_quadrature_at_edge[:,pt_index] = normal_vector
            end
        end

        dt_area = 0.0
        dt_centroid = [0.0,0.0]
        for pt_index in 1:length(quadrature_weight_at_edge)
            dt_area += 
                dot(u_g_quadrature_at_edge[:,pt_index], norm_quadrature_at_edge[:,pt_index]) * 
                    quadrature_weight_at_edge[:,pt_index]
            dt_centroid +=
                dot(u_g_quadrature_at_edge[:,pt_index], norm_quadrature_at_edge[:,pt_index]) * 
                    quadrature_weight_at_edge[:,pt_index] *
                    (geometry.centroid_coord - quadrature_coords_at_edge[:,pt_index])
        end
        dt_x_c = dt_centroid[1]
        dt_y_c = dt_centroid[2]



        return basis_dt_mat
    end

end # module MeshMotion
