module Meshes

using LinearAlgebra: norm


struct Face
    start_vertex_index::Int
    end_vertex_index::Int

    # Looking from start to end, the left cell lies on the left of the edge.
    # The right cell lies on the right; index 0 denotes an exterior boundary.
    left_cell_index::Int
    right_cell_index::Int

    contact::Bool

    boundary::Symbol
    periodic_shift_vec::Vector{Float64}
end

struct Mesh
    # Initial logical bounds: [x0, x1, y0, y1], or [r0, r1, theta0, theta1].
    initial_boundary_box::Vector{Float64}
    nx::Int # Number of cells along the first logical coordinate.
    ny::Int # Number of cells along the second logical coordinate.
    cells_vertex_indices_mat::Matrix{Int} # Four counterclockwise vertices per column.
    faces_list::Vector{Face}

    periodic::Bool
    vertex2group_index::Vector{Int} # Vertex index -> shared velocity-group index.

    polar::Bool
    has_origin::Bool
end

function vertex_number(i, j, nx)
    return i + 1 + (nx + 1) * j
end

function cell_number(i, j, nx)
    return i + nx * (j - 1)
end

function make_mesh(box, nx, ny;
                   polar = false,
                   periodic = false,
                   boundaries = (:wall, :wall, :wall, :wall),
                   saltzman = false)

    x0, x1, y0, y1 = box
    if saltzman
        @assert 0.0 < y1 - y0 < (x1 - x0) / pi "Saltzman initialization requires 0 < Ly < Lx/pi to guarantee a non-folding initial mesh."
    end

    number_of_vertices = (nx + 1) * (ny + 1)

    # Global mesh data.
    vertices_coordinates_mat = zeros(2, number_of_vertices)
    cells_vertex_indices_mat = zeros(Int, 4, nx * ny)
    vertex2group_index = zeros(Int, number_of_vertices)

    for j in 0:ny
        for i in 0:nx
            vertex_index = vertex_number(i, j, nx)
            logical_x = x0 + (x1 - x0) * i / nx
            logical_y = y0 + (y1 - y0) * j / ny

            if polar
                vertices_coordinates_mat[:, vertex_index] =
                    [logical_x * cos(logical_y), logical_x * sin(logical_y)]
            elseif saltzman
                vertices_coordinates_mat[:, vertex_index] =
                    [logical_x + (y1 - logical_y) * sin(pi * i / nx), logical_y]
            else
                vertices_coordinates_mat[:, vertex_index] = [logical_x, logical_y]
            end

            vertex2group_index[vertex_index] = vertex_index
            if periodic
                vertex2group_index[vertex_index] = 1 + mod(i, nx) + nx * mod(j, ny)
            elseif polar && x0 == 0.0 && i == 0
                # All logical copies of the polar origin share one velocity group.
                vertex2group_index[vertex_index] = 1
            end
        end
    end

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
            cell_vertices = cells_vertex_indices_mat[:, cell_index]

            # Right face of this logical cell.
            if i < nx
                push!(faces_list,
                      Face(cell_vertices[2], cell_vertices[3],
                                cell_index, cell_number(i + 1, j, nx),
                                false, :interior, [0.0, 0.0]))
            elseif periodic
                push!(faces_list,
                      Face(cell_vertices[2], cell_vertices[3],
                                cell_index, cell_number(1, j, nx),
                                false, :periodic, [-(x1 - x0), 0.0]))
            else
                push!(faces_list,
                      Face(cell_vertices[2], cell_vertices[3],
                                cell_index, 0,
                                false, boundaries[2], [0.0, 0.0]))
            end

            # Top face.
            if j < ny
                push!(faces_list,
                      Face(cell_vertices[3], cell_vertices[4],
                                cell_index, cell_number(i, j + 1, nx),
                                false, :interior, [0.0, 0.0]))
            elseif periodic
                push!(faces_list,
                      Face(cell_vertices[3], cell_vertices[4],
                                cell_index, cell_number(i, 1, nx),
                                false, :periodic, [0.0, -(y1 - y0)]))
            else
                push!(faces_list,
                      Face(cell_vertices[3], cell_vertices[4],
                                cell_index, 0,
                                false, boundaries[4], [0.0, 0.0]))
            end

            # Left and bottom boundary faces are added only once.
            if i == 1 && !periodic
                push!(faces_list,
                      Face(cell_vertices[4], cell_vertices[1],
                                cell_index, 0,
                                false, boundaries[1], [0.0, 0.0]))
            end
            if j == 1 && !periodic
                push!(faces_list,
                      Face(cell_vertices[1], cell_vertices[2],
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
        periodic,
        vertex2group_index,
        polar,
        polar && x0 == 0.0,
    )

    return mesh, vertices_coordinates_mat
end

function mark_contact_faces!(mesh::Mesh, vertices_coordinates_mat, is_contact_face)
    # Tag existing faces in physical coordinates, independently of their direction.
    # The tags stay with the connectivity while the mesh moves.
    for face_index in eachindex(mesh.faces_list)
        face = mesh.faces_list[face_index]
        contact = false
        if face.right_cell_index > 0
            start_point = vertices_coordinates_mat[:, face.start_vertex_index]
            end_point = vertices_coordinates_mat[:, face.end_vertex_index]
            contact = is_contact_face(start_point, end_point)
        end

        if contact != face.contact
            # I though it just C++ I can modify the attribute of a struct,
            # but in Julia, structs are immutable!!!!
            # . So I have to replace the whole struct in the vector!!!
            # Face is immutable;
            # replace only this vector entry to update its tag.
            mesh.faces_list[face_index] = Face(
                face.start_vertex_index, face.end_vertex_index,
                face.left_cell_index, face.right_cell_index,
                contact, face.boundary, face.periodic_shift_vec,
            )
        end
    end
    return mesh
end

function current_face_geometry(vertices_coordinates_mat, face::Face)
    start_point = vertices_coordinates_mat[:, face.start_vertex_index]
    end_point = vertices_coordinates_mat[:, face.end_vertex_index]

    edge_vector = end_point - start_point
    edge_length = norm(edge_vector)

    if edge_length < 1.0e-14
        return start_point,
               end_point,
               0.0,
               [0.0, 0.0]
    end

    # Rotate the directed edge vector clockwise by 90 degrees.
    # With the face orientation convention used here, this gives
    # the unit normal pointing from the left cell to the right cell.

    normal = [edge_vector[2], -edge_vector[1]] / edge_length

    return start_point,
           end_point,
           edge_length,
           normal
end

number_of_cells(mesh::Mesh) = size(mesh.cells_vertex_indices_mat, 2)
cell_vertex_indices(mesh::Mesh, cell_index::Int) = mesh.cells_vertex_indices_mat[:, cell_index]

end # module Meshes
