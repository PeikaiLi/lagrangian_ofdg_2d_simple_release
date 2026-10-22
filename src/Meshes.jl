module Meshes

    using LinearAlgebra: norm

    struct Face
        start_vertex_index::Int
        end_vertex_index::Int
        left_cell_index::Int
        right_cell_index::Int
        contact::Bool
        boundary::Symbol
        periodic_shift_vec::Vector{Float64}
    end

    struct Mesh
        initial_boundary_box::Vector{Float64}
        nx::Int 
        ny::Int 
        cells_vertex_indices_mat::Matrix{Int} 
        faces_list::Vector{Face}
        boundaries::NTuple{4,Symbol} # Left, right, bottom, top 
        vertex2group_index::Vector{Int} 
    end

    # I learned this mapping trick in 2025Fall from Per-olof, which is friendly for C++
    function vertex_number(i, j, nx)
        return i + 1 + (nx + 1) * j
    end

    function cell_number(i, j, nx)
        return i + nx * (j - 1)
    end

    function make_mesh(box, nx, ny; boundaries = (:wall, :wall, :wall, :wall))
        x0, x1, y0, y1 = box
        length(boundaries) == 4 || throw(ArgumentError("boundaries must contain left, right, bottom, and top"))
        boundaries = Tuple(boundaries)
        for boundary in boundaries
            boundary in (:wall, :fixed, :periodic, :outflow) ||
                throw(ArgumentError("unsupported boundary: $boundary"))
        end
        periodic_x = boundaries[1] == :periodic
        periodic_y = boundaries[3] == :periodic
        if periodic_x != (boundaries[2] == :periodic)
            throw(ArgumentError("left and right boundaries must both be periodic, or neither"))
        end
        if periodic_y != (boundaries[4] == :periodic)
            throw(ArgumentError("bottom and top boundaries must both be periodic, or neither"))
        end
        number_of_vertices = (nx + 1) * (ny + 1)
        number_of_group_columns = periodic_x ? nx : nx + 1

        # Global mesh data.
        vertices_coords_mat = zeros(2, number_of_vertices)
        cells_vertex_indices_mat = zeros(Int, 4, nx * ny)
        vertex2group_index = zeros(Int, number_of_vertices)

        for j in 0:ny
            for i in 0:nx
                vertex_index = vertex_number(i, j, nx)
                x_coord = x0 + (x1 - x0) * i / nx
                y_coord = y0 + (y1 - y0) * j / ny

                vertices_coords_mat[:, vertex_index] = [x_coord, y_coord]

                group_i = periodic_x ? mod(i, nx) : i
                group_j = periodic_y ? mod(j, ny) : j
                vertex2group_index[vertex_index] = 1 + group_i + number_of_group_columns * group_j
            end
        end

        # This index design is intended to have the relationship like this.
        for j in 1:ny
            for i in 1:nx
                cell_index = cell_number(i, j, nx)
                cells_vertex_indices_mat[:, cell_index] = [
                    vertex_number(i - 1, j - 1, nx),
                    vertex_number(i, j - 1, nx),
                    vertex_number(i, j, nx),
                    vertex_number(i - 1, j, nx),
                ]
            end
        end

        faces_list = Face[]
        for j in 1:ny
            for i in 1:nx
                cell_index = cell_number(i, j, nx)
                cell_vertex_indices_list = cells_vertex_indices_mat[:, cell_index]

                # Right face of this cell.
                if i < nx
                    push!(faces_list,
                        Face(cell_vertex_indices_list[2], cell_vertex_indices_list[3],
                                    cell_index, cell_number(i + 1, j, nx),
                                    false, :interior, [0.0, 0.0]))
                elseif periodic_x
                    push!(faces_list,
                        Face(cell_vertex_indices_list[2], cell_vertex_indices_list[3],
                                    cell_index, cell_number(1, j, nx),
                                    false, :periodic, [-(x1 - x0), 0.0]))
                else
                    push!(faces_list,
                        Face(cell_vertex_indices_list[2], cell_vertex_indices_list[3],
                                    cell_index, 0,
                                    false, boundaries[2], [0.0, 0.0]))
                end

                # Top face.
                if j < ny
                    push!(faces_list,
                        Face(cell_vertex_indices_list[3], cell_vertex_indices_list[4],
                                    cell_index, cell_number(i, j + 1, nx),
                                    false, :interior, [0.0, 0.0]))
                elseif periodic_y
                    push!(faces_list,
                        Face(cell_vertex_indices_list[3], cell_vertex_indices_list[4],
                                    cell_index, cell_number(i, 1, nx),
                                    false, :periodic, [0.0, -(y1 - y0)]))
                else
                    push!(faces_list,
                        Face(cell_vertex_indices_list[3], cell_vertex_indices_list[4],
                                    cell_index, 0,
                                    false, boundaries[4], [0.0, 0.0]))
                end

                # Left and bottom boundary faces are added only once.
                if i == 1 && !periodic_x
                    push!(faces_list,
                        Face(cell_vertex_indices_list[4], cell_vertex_indices_list[1],
                                    cell_index, 0,
                                    false, boundaries[1], [0.0, 0.0]))
                end
                if j == 1 && !periodic_y
                    push!(faces_list,
                        Face(cell_vertex_indices_list[1], cell_vertex_indices_list[2],
                                    cell_index, 0,
                                    false, boundaries[3], [0.0, 0.0]))
                end
            end
        end

        mesh = Mesh(
            collect(box),
            nx,
            ny,
            cells_vertex_indices_mat,
            faces_list,
            boundaries,
            vertex2group_index,
        )

        return mesh, vertices_coords_mat
    end

    function mark_contact_faces!(mesh::Mesh, vertices_coords_mat, is_contact_face)
        for face_index in eachindex(mesh.faces_list)
            face = mesh.faces_list[face_index]
            contact = false
            if face.right_cell_index > 0
                start_coords = vertices_coords_mat[:, face.start_vertex_index]
                end_coords = vertices_coords_mat[:, face.end_vertex_index]
                contact = is_contact_face(start_coords, end_coords)
            end

            if contact != face.contact
                # Face is immutable; Not like C++! Trickly, we can replace the face in the list with a new one.
                mesh.faces_list[face_index] = Face(
                    face.start_vertex_index, face.end_vertex_index,
                    face.left_cell_index, face.right_cell_index,
                    contact, face.boundary, face.periodic_shift_vec,
                )
            end
        end
        return mesh
    end

    function current_face_geometry(vertices_coords_mat, face::Face)
        # I do not want save this, since it will be changed in time.
        # If I rewrite it in C++, I will save it in the Face class.
        start_coords = vertices_coords_mat[:, face.start_vertex_index]
        end_coords = vertices_coords_mat[:, face.end_vertex_index]

        # from start to end
        edge_vector = end_coords - start_coords
        edge_length = norm(edge_vector)

        if edge_length  == 0.0
            return start_coords,
                end_coords,
                0.0,
                [0.0, 0.0]
        end

        # 90-degree clockwise rotation
        normal = [edge_vector[2], -edge_vector[1]] / edge_length

        return start_coords,
            end_coords,
            edge_length,
            normal
    end

    number_of_cells(mesh::Mesh) = size(mesh.cells_vertex_indices_mat, 2)
    cell_vertex_indices(mesh::Mesh, cell_index::Int) = mesh.cells_vertex_indices_mat[:, cell_index]

end # module Meshes
