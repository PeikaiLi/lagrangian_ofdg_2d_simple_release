module Quadrature

using LinearAlgebra: SymTridiagonal, eigen

# Refer to UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>

function gauss_legendre(number_of_points)
    if number_of_points == 1
        return [0.0], [2.0]
    end

    off_diagonal_list = zeros(number_of_points - 1)
    for point_index in 1:number_of_points - 1
        off_diagonal_list[point_index] = point_index / sqrt(4.0 * point_index^2 - 1.0)
    end

    jacobi_mat = SymTridiagonal(zeros(number_of_points), off_diagonal_list)
    eigen_decomposition = eigen(jacobi_mat)

    quadrature_points_1d_list = collect(eigen_decomposition.values)
    quadrature_weights_1d_list = 2.0 .* eigen_decomposition.vectors[1, :].^2

    return quadrature_points_1d_list, quadrature_weights_1d_list
end

function triangle_area_twice(vertex_a, vertex_b, vertex_c)
    return (vertex_b[1] - vertex_a[1]) * (vertex_c[2] - vertex_a[2]) -
           (vertex_b[2] - vertex_a[2]) * (vertex_c[1] - vertex_a[1])
end

function quadrilateral_quadrature(cell_vertices_mat,
                                  quadrature_points_1d_list,
                                  quadrature_weights_1d_list)

    first_twice_area = triangle_area_twice(
        cell_vertices_mat[:, 1],
        cell_vertices_mat[:, 2],
        cell_vertices_mat[:, 3],
    )
    second_twice_area = triangle_area_twice(
        cell_vertices_mat[:, 1],
        cell_vertices_mat[:, 3],
        cell_vertices_mat[:, 4],
    )

    x_span = maximum(cell_vertices_mat[1, :]) - minimum(cell_vertices_mat[1, :])
    y_span = maximum(cell_vertices_mat[2, :]) - minimum(cell_vertices_mat[2, :])
    tolerance = 1.0e-13 * max(x_span, y_span)^2

    if first_twice_area >= -tolerance && second_twice_area >= -tolerance
        two_triangle_vertex_indices_list = [(1, 2, 3), (1, 3, 4)]
    else
        first_twice_area = triangle_area_twice(
            cell_vertices_mat[:, 1],
            cell_vertices_mat[:, 2],
            cell_vertices_mat[:, 4],
        )
        second_twice_area = triangle_area_twice(
            cell_vertices_mat[:, 2],
            cell_vertices_mat[:, 3],
            cell_vertices_mat[:, 4],
        )
        if first_twice_area < -tolerance || second_twice_area < -tolerance
            throw(DomainError(first_twice_area + second_twice_area,
                              "crossed or inverted quadrilateral"))
        end
        two_triangle_vertex_indices_list = [(1, 2, 4), (2, 3, 4)]
    end

    number_of_1d_points = length(quadrature_points_1d_list)
    quadrature_point_vectors_list = Vector{Float64}[]
    quadrature_weights_list = Float64[]

    for triangle_vertices in two_triangle_vertex_indices_list
        vertex_a = cell_vertices_mat[:, triangle_vertices[1]]
        vertex_b = cell_vertices_mat[:, triangle_vertices[2]]
        vertex_c = cell_vertices_mat[:, triangle_vertices[3]]

        twice_area = triangle_area_twice(vertex_a, vertex_b, vertex_c)

        if twice_area <= 0.0
            continue # zero-area triangle at the polar origin
        end

        for first_quadrature_index in 1:number_of_1d_points
            duffy_u = 0.5 * (1.0 + quadrature_points_1d_list[first_quadrature_index])
            for second_quadrature_index in 1:number_of_1d_points
                duffy_v = 0.5 * (1.0 + quadrature_points_1d_list[second_quadrature_index])

                quadrature_point = vertex_a +
                    duffy_u * (vertex_b - vertex_a) + (1.0 - duffy_u) * duffy_v * (vertex_c - vertex_a)

                quadrature_weight = 0.25 * quadrature_weights_1d_list[first_quadrature_index] *
                    quadrature_weights_1d_list[second_quadrature_index] * (1.0 - duffy_u) * twice_area

                push!(quadrature_point_vectors_list, quadrature_point)
                push!(quadrature_weights_list, quadrature_weight)
            end
        end
    end

    if isempty(quadrature_weights_list)
        throw(DomainError(0.0, "zero cell area"))
    end

    quadrature_points_mat = hcat(quadrature_point_vectors_list...)
    return quadrature_points_mat, quadrature_weights_list
end

end # module Quadrature
