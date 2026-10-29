module FaceFlux

    using LinearAlgebra: I, dot

    import ..Euler
    import ..Geometry
    import ..Problems

    function outside_cons(inside_cons_vec, point_coord, boundary, normal, problem)
        if boundary == :wall
            return Euler.reflected_cons(inside_cons_vec, normal)
        elseif boundary == :fixed
            return Problems.initial_cons(problem, point_coord)
        elseif boundary == :outflow
            return copy(inside_cons_vec)
        end
        error("Unknown boundary: $(boundary)")
    end

    function face_cons(face, point_coord, cells_coefficients, 
        cells_geometries,
        problem::Problems.EulerProblem, 
        cells_gamma_ls, 
        exponents_ls, 
        normal
        )

        left_cell_index = face.left_cell_index
        right_cell_index = face.right_cell_index

        left_cons_vec = Geometry.eval_cons(
            cells_coefficients[:, :, left_cell_index],
            cells_geometries[left_cell_index],
            point_coord,
            exponents_ls,
        )
        left_gamma = cells_gamma_ls[left_cell_index]

        if right_cell_index > 0
            right_point_coord = point_coord + face.periodic_shift_vec
            right_cons_vec = Geometry.eval_cons(
                cells_coefficients[:, :, right_cell_index],
                cells_geometries[right_cell_index],
                right_point_coord,
                exponents_ls,
            )
            right_gamma = cells_gamma_ls[right_cell_index]
        else
            right_cons_vec = outside_cons(left_cons_vec, point_coord, face.boundary, normal, problem)
            right_gamma = left_gamma
        end

        return left_cons_vec,
            right_cons_vec,
            left_gamma,
            right_gamma
    end

    # Split a pointwise moving-face flux before any integration or DG mode selection.
    function moving_face_flux(left_cons_vec, right_cons_vec, left_gamma, right_gamma,
                            normal, mesh_velocity)
        star_normal_velocity, star_pressure = Euler.hllc_star(
            left_cons_vec, right_cons_vec, left_gamma, right_gamma, normal)
        relative_normal_speed = star_normal_velocity - dot(mesh_velocity, normal)

        left_pressure = Euler.pressure(left_cons_vec, left_gamma)
        right_pressure = Euler.pressure(right_cons_vec, right_gamma)
        if relative_normal_speed >= 0.0
            upwind_cons_vec = left_cons_vec
            upwind_pressure = left_pressure
            oriented_pressure_jump = left_pressure - right_pressure
        else
            upwind_cons_vec = right_cons_vec
            upwind_pressure = right_pressure
            oriented_pressure_jump = right_pressure - left_pressure
        end

        advection_flux = relative_normal_speed * upwind_cons_vec
        pressure_flux = [0.0, star_pressure * normal[1],
                        star_pressure * normal[2], star_pressure * star_normal_velocity]
        return advection_flux, pressure_flux, oriented_pressure_jump, upwind_pressure
    end

    # Extend boundary derivatives for the OF jump sensor. Coords are reflected
    # as well as momentum components at a stationary wall.
    function outside_cons_derivative(face,
                                point_coord,
                                normal,
                                cells_coefficients,
                                cells_geometries,
                                problem,
                                exponents_ls,
                                dx,
                                dy)

        cell_index = face.left_cell_index

        if dx + dy == 0
            inside_cons_vec = Geometry.eval_cons(
                cells_coefficients[:, :, cell_index],
                cells_geometries[cell_index],
                point_coord,
                exponents_ls,
            )
            return outside_cons(inside_cons_vec, point_coord, face.boundary, normal, problem)
        elseif face.boundary == :fixed
            return zeros(4)
        elseif face.boundary == :outflow
            return Geometry.eval_cons(
                cells_coefficients[:, :, cell_index],
                cells_geometries[cell_index],
                point_coord,
                exponents_ls;
                dx = dx,
                dy = dy,
            )
        end

        if face.boundary != :wall
            error("Unknown boundary: $(face.boundary)")
        end

        reflection_mat = Matrix{Float64}(I, 2, 2) - 2.0 * normal * transpose(normal)

        cons_derivative_vec = zeros(4)
        for x_derivative_index in 0:dx
            for y_derivative_index in 0:dy
                coefficient = binomial(dx, x_derivative_index) *
                    binomial(dy, y_derivative_index) * reflection_mat[1, 1]^x_derivative_index *
                    reflection_mat[2, 1]^(dx - x_derivative_index) *
                    reflection_mat[1, 2]^y_derivative_index * reflection_mat[2, 2]^(dy - y_derivative_index)

                inside_cons_derivative_vec = Geometry.eval_cons(
                    cells_coefficients[:, :, cell_index],
                    cells_geometries[cell_index],
                    point_coord,
                    exponents_ls;
                    dx = x_derivative_index + y_derivative_index,
                    dy = dx + dy - x_derivative_index - y_derivative_index,
                )
                cons_derivative_vec += coefficient * inside_cons_derivative_vec
            end
        end

        return Euler.reflected_cons(cons_derivative_vec, normal)
    end

end # module FaceFlux
