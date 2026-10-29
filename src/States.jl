module States

    import ..Basis
    import ..Geometry
    import ..Meshes
    import ..Options

    struct DGState
        # dim × number_of_vertices. We are 2d dim = 2; if I write 3D version, it is also fine
        vertices_coords_mat::Matrix{Float64}

        # number_of_basis × 4 × number_of_cells.
        cells_moments_array::Array{Float64, 3}
    end

    number_of_vertices(state::DGState) = size(state.vertices_coords_mat, 2)
    number_of_basis(state::DGState) = size(state.cells_moments_array, 1)

    cell_vertices_coords(state::DGState, mesh::Meshes.Mesh, cell_index::Int) =
        state.vertices_coords_mat[:, Meshes.cell_vertex_indices(mesh, cell_index)]

    function initial_State(vertices_coords_mat, mesh::Meshes.Mesh, options::Options.DGOptions, init_cons)
        cells_geometries = Geometry.build_cells_geometries(vertices_coords_mat, mesh, options)

        number_of_basis = length(Basis.basis_powers(options.degree))
        cells_moments_array = zeros(number_of_basis, 4, length(cells_geometries))

        for cell_index in eachindex(cells_geometries)
            geometry = cells_geometries[cell_index]

            for quadrature_point_index in eachindex(geometry.quadrature_weights_ls)
                quadrature_point_coord = geometry.quadrature_coords_ls[quadrature_point_index]
                cons_vec = init_cons(quadrature_point_coord, cell_index)

                for basis_index in 1:number_of_basis
                    cells_moments_array[basis_index, :, cell_index] +=
                        geometry.quadrature_weights_ls[quadrature_point_index] *
                        geometry.basis_at_quadrature_mat[basis_index, quadrature_point_index] * cons_vec
                end
            end
        end

        return DGState(copy(vertices_coords_mat), cells_moments_array)
    end

end # States