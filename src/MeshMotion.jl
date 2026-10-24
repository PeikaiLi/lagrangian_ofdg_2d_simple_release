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


    struct CellMotionCache
        cell_velocities_mat::Matrix{Float64} # same order and size as cell_vertices_coords_mat

        # Each column corresponds to one edge quadrature point.
        edge_quadrature_coords_mat::Matrix{Float64}

        edge_quadrature_weights_ls::Vector{Float64}
        edge_quadrature_velocities_mat::Matrix{Float64}
        edge_quadrature_normals_mat::Matrix{Float64}

        dt_area::Float64
        dt_x_c::Float64
        dt_y_c::Float64
        dt_hx::Float64
        dt_hy::Float64
    end


    function build_CellMotionCache(
        geometry::Geometry.CellGeometry,
        cell_velocities_mat::Matrix{Float64}, # same order and size as cell_vertices_coords_mat
        options::Options.DGOptions,
    )

        num_of_edge = size(geometry.cell_vertices_coords_mat)[2]
        quadrature1d_pts_ls,quadrature1d_wts_ls =  options.quadrature1d_pts_ls, options.quadrature1d_wts_ls

        edge_quadrature_coords_mat = zeros(2,num_of_edge*length(quadrature1d_pts_ls))
        edge_quadrature_weights_ls = zeros(num_of_edge*length(quadrature1d_pts_ls))
        edge_quadrature_velocities_mat = zeros(2,num_of_edge*length(quadrature1d_pts_ls))
        edge_quadrature_normals_mat = zeros(2,num_of_edge*length(quadrature1d_pts_ls))

        pt_index = 0
        for start_vex_index in 1:num_of_edge
            start_vex_coord = geometry.cell_vertices_coords_mat[:,start_vex_index]
            end_vex_index = mod1(start_vex_index+1,num_of_edge)
            end_vex_coord = geometry.cell_vertices_coords_mat[:,end_vex_index]

            start_vex_volcity = cell_velocities_mat[:,start_vex_index]
            end_vex_volcity = cell_velocities_mat[:,end_vex_index]

            tangent_vector = end_vex_coord - start_vex_coord
            length_tangent_vector = norm(tangent_vector)
            Options.check_positive(length_tangent_vector, "zero edge in moving basis")
            tangent_vector = tangent_vector/norm(tangent_vector) 
            normal_vector = [tangent_vector[2],-tangent_vector[1]] # clock 90 degree

            for i in 1:length(quadrature1d_pts_ls)
                pt_index += 1
                edge_quadrature_coords_mat[:,pt_index] = start_vex_coord + quadrature1d_pts_ls[i]*(end_vex_coord - start_vex_coord)
                edge_quadrature_weights_ls[pt_index] = quadrature1d_wts_ls[i] * length_tangent_vector
                edge_quadrature_velocities_mat[:,pt_index] = start_vex_volcity + quadrature1d_pts_ls[i]*(end_vex_volcity - start_vex_volcity)
                edge_quadrature_normals_mat[:,pt_index] = normal_vector
            end
        end

        dt_area = 0.0
        dt_centroid = [0.0,0.0]
        for pt_index in 1:length(edge_quadrature_weights_ls)
            dt_area += 
                dot(edge_quadrature_velocities_mat[:,pt_index], edge_quadrature_normals_mat[:,pt_index]) * 
                    edge_quadrature_weights_ls[pt_index]
            dt_centroid +=
                dot(edge_quadrature_velocities_mat[:,pt_index], edge_quadrature_normals_mat[:,pt_index]) * 
                    edge_quadrature_weights_ls[pt_index] *
                    (edge_quadrature_coords_mat[:,pt_index] - geometry.centroid_coord)
        end
        dt_centroid = dt_centroid/geometry.area
        dt_x_c = dt_centroid[1]
        dt_y_c = dt_centroid[2]

        dt_hx = cell_velocities_mat[1,geometry.xmax_vertex_index] - cell_velocities_mat[1,geometry.xmin_vertex_index]
        dt_hy = cell_velocities_mat[2,geometry.ymax_vertex_index] - cell_velocities_mat[2,geometry.ymin_vertex_index]        

        return CellMotionCache(
            cell_velocities_mat,
            edge_quadrature_coords_mat,
            edge_quadrature_weights_ls,
            edge_quadrature_velocities_mat,
            edge_quadrature_normals_mat,
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
        motioncache::CellMotionCache,
        exponents_ls::Vector{Tuple{Int,Int}},
        options::Options.DGOptions,
        pt_coord::Vector{Float64}
    )
               
        xi  = 2 * (pt_coord[1] - geometry.centroid_coord[1]) / geometry.hx
        eta = 2 * (pt_coord[2] - geometry.centroid_coord[2]) / geometry.hy
        
        dt_xi  = -2 * motioncache.dt_x_c/geometry.hx - 2 * (pt_coord[1] - geometry.centroid_coord[1])/(geometry.hx^2) * motioncache.dt_hx
        dt_eta = -2 * motioncache.dt_y_c/geometry.hy - 2 * (pt_coord[2] - geometry.centroid_coord[2])/(geometry.hy^2) * motioncache.dt_hy

        dt_raw_basis_vals = zeros(length(exponents_ls))
        for i in 1:length(exponents_ls)
            (a, b) = exponents_ls[i]
            k_l = Basis.monomial_factor(a, b)
            if a > 0 
                dt_raw_basis_vals[i] += 
                    a*(xi^(a-1)*eta^b)*dt_xi
            end
            if b > 0
                dt_raw_basis_vals[i] +=
                    b*(xi^a*eta^(b-1))*dt_eta
            end
            dt_raw_basis_vals[i]  = k_l * dt_raw_basis_vals[i]
        end
        return dt_raw_basis_vals

    end


    function raw_basis_time_derivative_at_quadratures(
        geometry::Geometry.CellGeometry,
        motioncache::CellMotionCache,
        exponents_ls::Vector{Tuple{Int,Int}},
        options::Options.DGOptions
    )

        dt_raw_basis_mat = zeros(size(geometry.basis_at_quadrature_mat))
        for i in 1:length(geometry.quadrature_coords_ls)
            pt_coord = geometry.quadrature_coords_ls[i]
            dt_raw_basis_mat[:, i] = raw_basis_time_derivative_at_one_pt(
                geometry,
                motioncache, 
                exponents_ls,
                options,
                pt_coord
            ) 
        end

        return dt_raw_basis_mat
    end


    function basis_time_derivative_at_one_pt(
        geometry::Geometry.CellGeometry,
        motioncache::CellMotionCache,
        exponents_ls::Vector{Tuple{Int,Int}},
        options::Options.DGOptions,
        pt_coord::Vector{Float64},
    )
    
        dt_raw_basis_mat = raw_basis_time_derivative_at_quadratures(
            geometry,
            motioncache, 
            exponents_ls,
            options
        ) 

        dt_raw_basis_means_ls = dt_raw_basis_mat * geometry.quadrature_weights_ls / geometry.area

        boundary_flux_raw_basis_means_ls = zeros(size(dt_raw_basis_means_ls))     
        for pt_index in 1:length(motioncache.edge_quadrature_weights_ls)
            edge_pt_coord = motioncache.edge_quadrature_coords_mat[:,pt_index]
            raw_basis_vals = Basis.eval_raw_basis(
                edge_pt_coord,
                geometry.centroid_coord,
                geometry.hx,
                geometry.hy,
                exponents_ls;
            )

            boundary_flux_raw_basis_means_ls +=
                raw_basis_vals * 
                motioncache.edge_quadrature_weights_ls[:,pt_index] *
                dot(motioncache.edge_quadrature_velocities_mat[:,pt_index],motioncache.edge_quadrature_normals_mat[:,pt_index])            
        end

        boundary_flux_raw_basis_means_ls = boundary_flux_raw_basis_means_ls/geometry.area

        dt_basis_vals =
            dt_raw_basis_means_ls + boundary_flux_raw_basis_means_ls - geometry.raw_basis_means_ls/geometry.area * motioncache.dt_area


        dt_raw_basis_vals = raw_basis_time_derivative_at_one_pt(
                geometry,
                motioncache, 
                exponents_ls,
                options,
                pt_coord
            ) 
        
        dt_basis_vals = dt_raw_basis_vals - dt_basis_vals

        return dt_basis_vals
    end 

    function basis_time_derivative_at_quadratures(
        geometry::Geometry.CellGeometry,
        motioncache::CellMotionCache,
        exponents_ls::Vector{Tuple{Int,Int}},
        options::Options.DGOptions
    )

        dt_basis_mat = zeros(size(geometry.basis_at_quadrature_mat))
        for i in 1:length(geometry.quadrature_coords_ls)
            pt_coord = geometry.quadrature_coords_ls[i]
            dt_basis_mat[:, i] = basis_time_derivative_at_one_pt(
                geometry,
                motioncache, 
                exponents_ls,
                options,
                pt_coord
            ) 
        end

        return dt_basis_mat
    end


end # module MeshMotion
