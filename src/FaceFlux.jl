module FaceFlux

using LinearAlgebra: dot
import ..Euler
import ..Geometry
import ..Problems

function is_wall(boundary)
    return boundary in (:wall, :piston, :circular_wall)
end

function outside_state(inside_state, point_coordinates, time, face, normal, problem)
    if is_wall(face.boundary)
        wall_velocity_vec = Problems.wall_velocity(problem, face.boundary, point_coordinates, time)
        wall_normal_speed = dot(wall_velocity_vec, normal)
        return Euler.reflected_state(inside_state, normal, wall_normal_speed)
    elseif face.boundary == :fixed
        return Problems.initial_state(problem, point_coordinates, face.left_cell_index)
    elseif face.boundary in (:outflow, :pressure_inner, :pressure_outer, :origin)
        return copy(inside_state)
    end

    error("unknown boundary condition")
end


function face_states(face,
                     point_coordinates,
                     time,
                     cells_coefficients,
                     cells_geometries,
                     problem,
                     basis_power_pairs_list,
                     normal)

    left_cell_index = face.left_cell_index
    right_cell_index = face.right_cell_index

    left_state = Geometry.evaluate_solution(
        cells_coefficients[left_cell_index],
        cells_geometries[left_cell_index],
        point_coordinates,
        basis_power_pairs_list,
    )
    left_gamma = problem.cells_gamma_list[left_cell_index]

    if right_cell_index > 0
        right_point_coordinates = point_coordinates + face.periodic_shift_vec
        right_state = Geometry.evaluate_solution(
            cells_coefficients[right_cell_index],
            cells_geometries[right_cell_index],
            right_point_coordinates,
            basis_power_pairs_list,
        )
        right_gamma = problem.cells_gamma_list[right_cell_index]
    else
        right_state = outside_state(left_state, point_coordinates, time, face, normal, problem)
        right_gamma = left_gamma
    end

    return left_state,
           right_state,
           left_gamma,
           right_gamma
end

function face_star(left_state, right_state, left_gamma, right_gamma, face, normal, problem, time)

    if face.boundary in (:pressure_inner, :pressure_outer)
        rho, velocity_x, velocity_y, pressure_value = Euler.cons_to_prim(left_state, left_gamma)
        normal_velocity = velocity_x * normal[1] + velocity_y * normal[2]
        prescribed_pressure = Problems.boundary_pressure(problem, face.boundary, time)
        star_normal_velocity = normal_velocity - Euler.wave_curve(
                prescribed_pressure,
                rho,
                pressure_value,
                left_gamma,
            )
        return star_normal_velocity, prescribed_pressure
    end

    return Euler.hllc_star(left_state, right_state, left_gamma, right_gamma, normal)
end

end # module FaceFlux
