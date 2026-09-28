module States

import ..Basis
import ..Geometry
import ..Meshes

# The state stores coordinates and moments; projection accepts an initial-state function.
struct DGState
    # 2 × number_of_vertices. Column `vertex_index` is one physical vertex.
    vertices_coordinates_mat::Matrix{Float64}

    # number_of_basis_functions × 4 × number_of_cells.
    # These are integral moments, not point values or polynomial coefficients.
    cells_moments_array::Array{Float64, 3}
end

function project_initial_condition(vertices_coordinates_mat, mesh, options, initial_state)
    cells_geometries = Geometry.build_cells_geometries(vertices_coordinates_mat, mesh, options)

    number_of_basis_functions = length(Basis.basis_powers(options.degree))
    cells_moments_array = zeros(number_of_basis_functions, 4, length(cells_geometries))

    for cell_index in eachindex(cells_geometries)
        geometry = cells_geometries[cell_index]

        for quadrature_point_index in eachindex(geometry.quadrature_weights_list)
            quadrature_point_coordinates = geometry.quadrature_coordinates_mat[:, quadrature_point_index]
            state_vec = initial_state(quadrature_point_coordinates, cell_index)

            for basis_index in 1:number_of_basis_functions
                cells_moments_array[basis_index, :, cell_index] +=
                    geometry.quadrature_weights_list[quadrature_point_index] *
                    geometry.basis_at_quadrature_mat[basis_index, quadrature_point_index] * state_vec
            end
        end
    end

    return DGState(copy(vertices_coordinates_mat), cells_moments_array)
end

number_of_vertices(state::DGState) = size(state.vertices_coordinates_mat, 2)
number_of_basis_functions(state::DGState) = size(state.cells_moments_array, 1)
cell_vertex_coordinates(state::DGState, mesh::Meshes.Mesh, cell_index::Int) =
    state.vertices_coordinates_mat[:, Meshes.cell_vertex_indices(mesh, cell_index)]

end # module States
