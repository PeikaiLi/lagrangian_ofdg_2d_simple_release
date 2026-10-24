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


    struct Cell_move_util
        quadrature_coords_at_edge,
        quadrature_weight_at_edge,
        u_g_quadrature_at_edge,
        norm_quadrature_at_edge,
        dt_area,
        dt_x_c,
        dt_y_c,
        dt_hx,
        dt_hy
    end


    function build_Cell_move_util(
        geometry::Geometry.CellGeometry,
        cell_velocities_mat::Matrix{Float64}, # same order and size as cell_vertices_coords_mat
        options::Options.DGOptions,
    )

        num_of_edge = size(geometry.cell_vertices_coords_mat)[2]
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

        dt_hx = geometry.extreme_vertex_indices_ls[2] - geometry.extreme_vertex_indices_ls[1] 
        dt_hy = geometry.extreme_vertex_indices_ls[4] - geometry.extreme_vertex_indices_ls[3]
        # extreme_vertex_indices_ls is
        # [
        #     xmin_vertex_index,
        #     xmax_vertex_index,
        #     ymin_vertex_index,
        #     ymax_vertex_index,
        # ]

        return Cell_move_util(
            quadrature_coords_at_edge,
            quadrature_weight_at_edge,
            u_g_quadrature_at_edge,
            norm_quadrature_at_edge,
            dt_area,
            dt_x_c,
            dt_y_c,
            dt_hx,
            dt_hy
        )

    end


    # raw_basis_time_derivative
    function raw_basis_time_derivative_at_one_pt(
        geometry::Geometry.CellGeometry,
        cell_velocities_mat::Matrix{Float64}, # same order and size as cell_vertices_coords_mat
        exponents_ls::Vector{Tuple{Int,Int}},
        options::Options.DGOptions,
        pt_coord::Vector{Float64}
    )
               
        basis_vals = eval_basis_values(geometry, pt_coord, exponents_ls) 
        raw_basis_vals = basis_vals + geometry.raw_basis_means_ls
        xi = raw_basis_vals[2]
        eta = raw_basis_vals[3]

        

             



        dt_xi = -2 * dt_x_c/geometry.hx - 2 * (pt_coord[1] - geometry.centroid_coord[1])/(geometry.hx^2) * dt_hx
        dt_eta = -2 * dt_y_c/geometry.hy - 2 * (pt_coord[2] - geometry.centroid_coord[2])/(geometry.hy^2) * dt_hy

        dt_raw_basis_vals = zeros(length(exponents_ls))
        for i in 1:length(exponents_ls)
            (a, b) = exponents_ls[i]
            k_l = monomial_factor(a, b)
            dt_raw_basis_vals[i] = 
                a*(xi^(a-1)*eta^b)*dt_xi +
                b*(xi^a*eta^(b-1))*dt_eta
            dt_raw_basis_vals[i]  = k_l * dt_raw_basis_vals[i]
        end
        return dt_raw_basis_vals
    end


    function raw_basis_time_derivative_at_quadratures(
        geometry::Geometry.CellGeometry,
        cell_velocities_mat::Matrix{Float64}, # same order and size as cell_vertices_coords_mat
        exponents_ls::Vector{Tuple{Int,Int}},
        options::Options.DGOptions
    )

        dt_raw_basis_mat = zeros(size(geometry.basis_at_quadrature_mat))
        for i in 1:length(geometry.quadrature_coords_ls)
            pt_coord = geometry.quadrature_coords_ls[i]
            dt_raw_basis_mat[:, i] = raw_basis_time_derivative_at_one_pt(
                geometry,
                cell_velocities_mat, 
                exponents_ls,
                options,
                pt_coord
            ) 
        end

        return dt_raw_basis_mat
    end


    function basis_time_derivative_at_one_pt(
        geometry::Geometry.CellGeometry,
        cell_velocities_mat::Matrix{Float64}, # same order and size as cell_vertices_coords_mat
        exponents_ls::Vector{Tuple{Int,Int}},
        options::Options.DGOptions
        pt_coord::Vector{Float64},
    )

        quadrature_coords_at_edge,
            quadrature_weight_at_edge,
            u_g_quadrature_at_edge,
            norm_quadrature_at_edge = boundary_integral_with_ug(
            geometry,
            cell_velocities_mat,
            options
            )

        dt_raw_basis_mat = raw_basis_time_derivative_at_quadratures(
            geometry,
            cell_velocities_mat, 
            exponents_ls,
            options
        ) 

        dt_raw_basis_means_ls = dt_raw_basis_mat * geometry.quadrature_weights_ls / geometry.area

        dt_basis_ls = dt_raw_basis_means_ls

        dt_basis_ls = dt_basis_ls +  



         
        for i in 1:length(quadrature_coords_ls)
            pt_coord = quadrature_coords_ls[i]
            dt_raw_basis_mat[:, i] = raw_basis_time_derivative_at_one_pt(
                geometry,
                cell_velocities_mat, 
                exponents_ls,
                pt_coord,
                options
            ) 
        end

        return dt_basis_mat
    end 


end # module MeshMotion
