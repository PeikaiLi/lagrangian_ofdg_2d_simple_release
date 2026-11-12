module OF

    import ..FaceFlux, ..Geometry, ..Meshes, ..States
    import ..Discretization
    import ..Basis, ..Euler


    function cells_jump_V(
        state::States.DGState, 
        case::Discretization.DGCase, 
        cells_geometries, 
        cells_coefficients,
        dx,dy)

        degree = case.options.degree
        exponents_ls = Basis.basis_powers(degree)
        number_of_cells = length(cells_geometries)

        jump_V = zeros(4,number_of_cells)

        for face in case.mesh.faces_ls

            start_coord, end_coord,
                edge_length,
                normal = Meshes.current_face_geometry(state.vertices_coords_mat, face)

            left = face.left_cell_index
            right = face.right_cell_index

            for pt in [start_coord, end_coord]
                DU_left = Geometry.eval_cons(cells_coefficients[:, :, left], cells_geometries[left], pt, exponents_ls; dx = dx, dy = dy)
                                
                left_ave_cons_vec = cells_coefficients[1, :, left]
                left_gamma = case.cells_gamma_ls[left]
                R_left = Euler.left_eigen_vector_Jacobians(left_ave_cons_vec, left_gamma, normal)


                if right > 0  
                    DU_right = Geometry.eval_cons(cells_coefficients[:, :, right], cells_geometries[right], pt + face.periodic_shift_vec, exponents_ls; dx = dx, dy = dy)
                    right_ave_cons_vec = cells_coefficients[1, :, right]
                    right_gamma = case.cells_gamma_ls[right]
                    R_right = Euler.left_eigen_vector_Jacobians(right_ave_cons_vec, right_gamma, -normal)
                    jump = (DU_left .- DU_right)
                    left_jump = R_left * jump
                    right_jump = R_right * jump

                    jump_V[:, left] += left_jump.^2
                    jump_V[:, right] += right_jump.^2

                else
                    DU_right = FaceFlux.outside_cons_derivative(face,
                        pt,
                        normal,
                        cells_coefficients,
                        cells_geometries,
                        case.problem,
                        exponents_ls,
                        dx,
                        dy)

                    jump = (DU_left .- DU_right)
                    left_jump = R_left * jump
                    jump_V[:, left] += left_jump.^2
                end   

            end
        end

        jump_V = jump_V ./ 4
        return jump_V.^(1/2)
    end



    function damping_coefficients(
            state::States.DGState, 
            case::Discretization.DGCase,                                
            cells_geometries, 
            cells_coefficients
        )

        degree = case.options.degree
        exponents_ls = Basis.basis_powers(degree)
        number_of_cells = length(cells_geometries)

        sigma_mat = zeros(degree + 1, number_of_cells)

        if degree == 0 
            return sigma_mat
        end

        h_tilde_ls = [geometry.longest_edge for geometry in cells_geometries]
        for face in case.mesh.faces_ls
            left = face.left_cell_index
            right = face.right_cell_index
            _, _, face_length, normal = Meshes.current_face_geometry(state.vertices_coords_mat, face)
            face_length > 0 || throw(DomainError(face_length, "OF requires nonzero edges"))

            if right > 0
                h_tilde_ls[left] = max(h_tilde_ls[left], cells_geometries[right].longest_edge)
                h_tilde_ls[right] = max(h_tilde_ls[right], cells_geometries[left].longest_edge)
            end
        end


        for a in 0:degree
            jump_mat = zeros(4,number_of_cells)
            derivatives_ls = [(dx, a - dx) for dx in 0:a]
            for (dx, dy) in derivatives_ls
                jump_mat .+= cells_jump_V(state, case, cells_geometries, cells_coefficients, dx, dy)
            end

            max_jump = vec(maximum(jump_mat; dims=1))
            max_jump = max_jump * (2.0 * (2*a + 1)) / ( (2*degree - 1) * factorial(a) )
            max_jump = max_jump .* h_tilde_ls.^(a)
            sigma_mat[a + 1, :] = max_jump
        end

        return sigma_mat
    end

    

end