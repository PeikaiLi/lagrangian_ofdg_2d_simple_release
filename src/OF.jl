module OF


    import ..Discretization
    import ..Basis, ..Euler, ..FaceFlux, ..Geometry, ..Meshes, ..States

    function damping_coefficients(state::States.DGState, case::Discretization.DGCase,
                                cells_geometries, cells_coefficients)

        degree = case.options.degree
        number_of_cells = length(cells_geometries)

        sigma_mat = zeros(degree + 1, number_of_cells)

        if degree == 0 
            return sigma_mat
        end

        # exponents_ls = Basis.basis_powers(degree)
        # derivatives_ls = [(dx, a - dx) for a in 0:degree for dx in 0:a]
        # jump_squares_array = zeros(4, length(derivatives_ls), number_of_cells)

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


    end



end