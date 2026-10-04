module MeshMotion

using LinearAlgebra: dot, norm
import ..Basis
import ..Euler
import ..Problems
import ..FaceFlux
import ..Meshes

function compute_vertex_velocities(state,
                                   mesh,
                                   problem,
                                   options,
                                   time,
                                   cells_geometries,
                                   cells_coefficients)

    basis_power_pairs_list = Basis.basis_powers(options.degree)
    number_of_velocity_groups = maximum(mesh.vertex2group_index)

    group_velocity_sums_mat = zeros(2, number_of_velocity_groups)
    group_weight_sums_list = zeros(number_of_velocity_groups)
    boundary_group_flags_list = falses(number_of_velocity_groups)

    for face in mesh.faces_list
        if face.right_cell_index == 0 && FaceFlux.is_wall(face.boundary)
            start_group_index = mesh.vertex2group_index[face.start_vertex_index]
            end_group_index = mesh.vertex2group_index[face.end_vertex_index]
            boundary_group_flags_list[start_group_index] = true
            boundary_group_flags_list[end_group_index] = true
        end
    end

    for face in mesh.faces_list
        _, _, face_length, normal = Meshes.current_face_geometry(state.vertices_coordinates_mat, face)
        if face_length == 0.0
            continue
        end

        tangent = [-normal[2], normal[1]]

        for vertex_index in (face.start_vertex_index, face.end_vertex_index)
            vertex_coordinates = state.vertices_coordinates_mat[:, vertex_index]

            left_state,
            right_state,
            left_gamma,
            right_gamma = FaceFlux.face_states(
                face,
                vertex_coordinates,
                time,
                cells_coefficients,
                cells_geometries,
                problem,
                basis_power_pairs_list,
                normal,
            )

            star_normal_velocity, _ = FaceFlux.face_star(
                left_state,
                right_state,
                left_gamma,
                right_gamma,
                face,
                normal,
                problem,
                time,
            )

            left_velocity = left_state[2:3] / left_state[1]
            right_velocity = right_state[2:3] / right_state[1]

            tangent_velocity = dot(0.5 * (left_velocity + right_velocity), tangent)

            candidate_velocity = star_normal_velocity * normal + tangent_velocity * tangent

            group_index = mesh.vertex2group_index[vertex_index]
            weight = 1.0
            if boundary_group_flags_list[group_index] && face.right_cell_index > 0
                weight = 2.0
            end

            group_velocity_sums_mat[:, group_index] += weight * candidate_velocity
            group_weight_sums_list[group_index] += weight
        end
    end

    vertices_velocities_mat = zeros(size(state.vertices_coordinates_mat))
    for vertex_index in 1:size(state.vertices_coordinates_mat, 2)
        group_index = mesh.vertex2group_index[vertex_index]
        if group_weight_sums_list[group_index] > 0.0
            vertices_velocities_mat[:, vertex_index] = group_velocity_sums_mat[:, group_index] /
                group_weight_sums_list[group_index]
        end
    end

    # Project wall vertices so their normal mesh speed matches the wall speed.
    for face in mesh.faces_list
        if face.right_cell_index == 0 && FaceFlux.is_wall(face.boundary)
            _, _, face_length, normal = Meshes.current_face_geometry(state.vertices_coordinates_mat, face)
            if face_length == 0.0
                continue
            end

            for vertex_index in (face.start_vertex_index, face.end_vertex_index)
                wall_normal = normal
                if face.boundary == :circular_wall
                    wall_normal = state.vertices_coordinates_mat[:, vertex_index] /
                        norm(state.vertices_coordinates_mat[:, vertex_index])
                end

                # wall_speed is defined depends on question case by case. SO 
                # I put this function in SpecialProblems.jl.
                # Only used in Saltzmann problem. Rest of the problems so far, wall speed is 0.0.
                target_normal_speed = dot(Problems.wall_velocity(problem, face.boundary, state.vertices_coordinates_mat[:, vertex_index], time), wall_normal)

                vertices_velocities_mat[:, vertex_index] += (target_normal_speed -
                     dot(vertices_velocities_mat[:, vertex_index], wall_normal)) * wall_normal
            end
        end
    end

    if mesh.has_origin
        for vertex_index in 1:size(state.vertices_coordinates_mat, 2)
            if mesh.vertex2group_index[vertex_index] == 1
                vertices_velocities_mat[:, vertex_index] .= 0.0
            end
        end
    end

    return vertices_velocities_mat
end

function basis_time_derivative(geometry,
                               cell_velocities_mat,
                               basis_power_pairs_list,
                               reference_quadrature_coordinates_list,
                               quadrature_weights_1d_list)

    area_rate = 0.0
    centroid_rate = zeros(2)
    boundary_integrals = zeros(length(basis_power_pairs_list))

    for local_edge_index in 1:4
        next_local_vertex_index = mod1(local_edge_index + 1, 4)

        edge_start_coordinates = geometry.vertex_coordinates_mat[:, local_edge_index]
        edge_end_coordinates = geometry.vertex_coordinates_mat[:, next_local_vertex_index]
        edge_vector = edge_end_coordinates - edge_start_coordinates
        edge_length = norm(edge_vector)

        if edge_length < 1.0e-14
            continue
        end

        edge_normal = [edge_vector[2], -edge_vector[1]] / edge_length

        for quadrature_point_index in eachindex(reference_quadrature_coordinates_list)
            edge_coordinate = 0.5 * (1.0 + reference_quadrature_coordinates_list[quadrature_point_index])

            quadrature_point_coordinates = (1.0 - edge_coordinate) * edge_start_coordinates + edge_coordinate * edge_end_coordinates

            grid_velocity = (1.0 - edge_coordinate) *
                cell_velocities_mat[:, local_edge_index] + edge_coordinate *
                cell_velocities_mat[:, next_local_vertex_index]

            boundary_weight = 0.5 * edge_length *
                quadrature_weights_1d_list[quadrature_point_index] * dot(grid_velocity, edge_normal)

            area_rate += boundary_weight
            centroid_rate += boundary_weight * (quadrature_point_coordinates - geometry.centroid_coordinates)
            boundary_integrals += boundary_weight * Basis.raw_basis(
                    quadrature_point_coordinates,
                    geometry.centroid_coordinates,
                    geometry.hx,
                    geometry.hy,
                    basis_power_pairs_list,
                )
        end
    end

    centroid_rate /= geometry.area

    xmin_vertex_index,
    xmax_vertex_index,
    ymin_vertex_index,
    ymax_vertex_index = geometry.scale_vertex_indices_list

    hx_rate = cell_velocities_mat[1, xmax_vertex_index] - cell_velocities_mat[1, xmin_vertex_index]
    hy_rate = cell_velocities_mat[2, ymax_vertex_index] - cell_velocities_mat[2, ymin_vertex_index]

    raw_basis_dt_mat = zeros(length(basis_power_pairs_list), length(geometry.quadrature_weights_list))

    for quadrature_point_index in eachindex(geometry.quadrature_weights_list)
        quadrature_point_coordinates = geometry.quadrature_coordinates_mat[:, quadrature_point_index]

        x_rate_at_fixed_point = -centroid_rate[1] -
            (quadrature_point_coordinates[1] - geometry.centroid_coordinates[1]) * hx_rate / geometry.hx

        y_rate_at_fixed_point = -centroid_rate[2] -
            (quadrature_point_coordinates[2] - geometry.centroid_coordinates[2]) * hy_rate / geometry.hy

        raw_basis_dt_mat[:, quadrature_point_index] = x_rate_at_fixed_point *
            geometry.basis_dx_at_quadrature_mat[:, quadrature_point_index] + y_rate_at_fixed_point *
            geometry.basis_dy_at_quadrature_mat[:, quadrature_point_index]
    end

    mean_rates = (raw_basis_dt_mat * geometry.quadrature_weights_list +
         boundary_integrals - geometry.raw_means_list * area_rate) / geometry.area

    basis_dt_mat = zeros(size(raw_basis_dt_mat))
    for quadrature_point_index in eachindex(geometry.quadrature_weights_list)
        basis_dt_mat[:, quadrature_point_index] = raw_basis_dt_mat[:, quadrature_point_index] -
            mean_rates
        basis_dt_mat[1, quadrature_point_index] = 0.0
    end

    return basis_dt_mat
end

end # module MeshMotion
