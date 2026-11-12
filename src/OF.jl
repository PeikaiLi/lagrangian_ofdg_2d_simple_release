module OF

    using LinearAlgebra

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

            if case.options.use_characteristic_decomposition                
                left_ave_cons_vec = cells_coefficients[1, :, left]
                left_gamma = case.cells_gamma_ls[left]
                R_left = Euler.left_eigen_vector_Jacobians(left_ave_cons_vec, left_gamma, normal)
                if right > 0
                    right_ave_cons_vec = cells_coefficients[1, :, right]
                    right_gamma = case.cells_gamma_ls[right]
                    R_right = Euler.left_eigen_vector_Jacobians(right_ave_cons_vec, right_gamma, -normal)
                end
            end

            for pt in [start_coord, end_coord]
                DU_left = Geometry.eval_cons(cells_coefficients[:, :, left], cells_geometries[left], pt, exponents_ls; dx = dx, dy = dy)

                if right > 0  
                    DU_right = Geometry.eval_cons(cells_coefficients[:, :, right], 
                        cells_geometries[right], 
                        pt + face.periodic_shift_vec, 
                        exponents_ls; 
                        dx = dx, dy = dy)
                    jump = (DU_left .- DU_right)
                    if case.options.use_characteristic_decomposition
                        right_jump = R_right * jump
                        jump_V[:, right] += right_jump.^2
                        left_jump = R_left * jump
                        jump_V[:, left] += left_jump.^2
                    else
                        jump_V[:, left] += jump.^2
                        jump_V[:, right] += jump.^2
                    end

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

                    if case.options.use_characteristic_decomposition
                        left_jump = R_left * jump
                        jump_V[:, left] += left_jump.^2
                    else
                        jump_V[:, left] += jump.^2
                    end
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


    function build_project_mat(cell_mass_mat, project_degree, exponents_ls)

        number_of_basis = length(exponents_ls)
        project_mat = zeros(Float64,number_of_basis,number_of_basis)

        if project_degree == -1
            project_degree = 0
        end

        indice = 0
        for power in exponents_ls
            if sum(power) <= project_degree
                indice += 1
            end
        end

        if indice == 0
            return project_mat
        end

        project_mat[1:indice, :] = 
            cell_mass_mat[1:indice, 1:indice] \ cell_mass_mat[1:indice, :]

        return project_mat
        
     end


    function build_damping_mat(
            geometry::Meshes.CellGeometry,
            sigma_ls, # sigma_mat[:, cell_index]
            exponents_ls)

        number_of_basis = length(exponents_ls)
        id_mat = Matrix{Float64}(I, number_of_basis, number_of_basis)

        h = geometry.longest_edge
        cell_mass_mat = geometry.mass_mat
        damping_mat = zeros(Float64,number_of_basis,number_of_basis)

        for a in eachindex(sigma_ls)
            damping_mat +=
                sigma_ls[a] * cell_mass_mat * (
                    id_mat - build_project_mat(cell_mass_mat, a - 2, exponents_ls)
                    )
        end
        
        damping_mat = damping_mat / h
        damping_mat = cell_mass_mat \ damping_mat
        
        return damping_mat
     end


    function apply_filter(state::States.DGState, case::Discretization.DGCase, dt::Real)

        (case.options.degree == 0 || dt == 0) && return state

        exponents_ls = Basis.basis_powers(case.options.degree)

        cells_geometries = Geometry.build_cells_geometries(state.vertices_coords_mat, case.mesh, case.options)
        coefficients = Geometry.moments2coefficients(state.cells_moments_array, cells_geometries)

        sigma_mat = damping_coefficients(state, case, cells_geometries, coefficients)
        
        filtered_coefficients = copy(coefficients)
        for cell_index in eachindex(cells_geometries)

            damping_mat = build_damping_mat(cells_geometries[cell_index], sigma_mat[:, cell_index], exponents_ls)

            # exp(matrix), not exp.(matrix)
            filtered_coefficients[:, :, cell_index] =
                exp(-dt * damping_mat) * coefficients[:, :, cell_index]

            # for test
            # c0_before = coefficients[1, :, cell_index]
            # c0_after = filtered_coefficients[1, :, cell_index]

            # if !all(isapprox.(c0_after, c0_before; rtol=1e-10, atol=1e-12))
            #     error("OF changed the cell average in cell $cell_index")
            # end

        end

        moments = Geometry.coefficients2moments(filtered_coefficients, cells_geometries)

        if !all(isapprox.(moments[1, :, :], state.cells_moments_array[1, :, :]; rtol=1e-10, atol=1e-12))
            error("OF changed the cell average")
        end

        moments[1, :, :] = state.cells_moments_array[1, :, :]

        return States.DGState(state.vertices_coords_mat, moments)
    end

    
end