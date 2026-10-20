module Discretization

    import ..Meshes
    import ..Options
    import ..Problems
    import ..States

    struct DGCase{P<:Problems.EulerProblem}

        problem::P
        mesh::Meshes.Mesh
        options::Options.DGOptions

        initial::States.DGState
        cells_gamma_list::Vector{Float64}
    end

    function build_cells_gamma(problem::Problems.EulerProblem, mesh::Meshes.Mesh, init_vex_coords_mat)
        cells_gamma_list = zeros(Meshes.number_of_cells(mesh))
        for cell_index in eachindex(cells_gamma_list)
            cell_vertices = mesh.cells_vertex_indices_mat[:, cell_index]
            center_coordinates = zeros(2)
            for vertex_index in cell_vertices
                center_coordinates += init_vex_coords_mat[:, vertex_index]
            end
            center_coordinates /= length(cell_vertices)
            gamma = Problems.initial_gamma(problem, center_coordinates)
            if !isfinite(gamma) || gamma <= 1.0
                throw(ArgumentError("initial material gamma must be finite and greater than one"))
            end
            cells_gamma_list[cell_index] = gamma
        end
        return cells_gamma_list
    end

    function DGCase(problem::Problems.EulerProblem, mesh::Meshes.Mesh,
                    init_vex_coords_mat, options::Options.DGOptions;
                    init_cons = (point_coords, cell_index) ->
                        Problems.initial_cons(problem, point_coords))
        cells_gamma_list = build_cells_gamma(problem, mesh, init_vex_coords_mat)
        initial = States.initial_State(init_vex_coords_mat, mesh, options, init_cons)
        return DGCase(problem, mesh, options, initial, cells_gamma_list)
    end

end # module Discretization