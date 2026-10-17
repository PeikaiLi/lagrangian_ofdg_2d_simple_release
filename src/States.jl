module States

import ..Basis
import ..Geometry
import ..Meshes

    struct DGState
        # dim × number_of_vertices. We are 2d dim = 2; if I write 3D version, it is also fine
        vertices_coordinates_mat::Matrix{Float64}

        # number_of_basis × 4 × number_of_cells.
        cells_moments_array::Array{Float64, 3}
    end

    number_of_vertices(state::DGState) = size(state.vertices_coordinates_mat, 2)
    number_of_basis(state::DGState) = size(state.cells_moments_array, 1)

    cell_vertex_coordinates(state::DGState, mesh::Meshes.Mesh, cell_index::Int) =
        state.vertices_coordinates_mat[:, Meshes.cell_vertex_indices(mesh, cell_index)]

end # States


