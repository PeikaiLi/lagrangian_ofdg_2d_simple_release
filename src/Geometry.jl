module Geometry

using LinearAlgebra: norm
import ..Basis
import ..Meshes
import ..Options
import ..Quadrature

# Current-cell quadrature, basis values, and the mass matrix are rebuilt after mesh motion.
struct CellGeometry
    vertex_coordinates_mat::Matrix{Float64}
    quadrature_coordinates_mat::Matrix{Float64}
    quadrature_weights_list::Vector{Float64}
    area::Float64
    centroid_coordinates::Vector{Float64}
    hx::Float64
    hy::Float64
    scale_vertex_indices_list::Vector{Int}
    raw_means_list::Vector{Float64}
    mass_mat::Matrix{Float64}
    basis_at_quadrature_mat::Matrix{Float64}
    basis_dx_at_quadrature_mat::Matrix{Float64}
    basis_dy_at_quadrature_mat::Matrix{Float64}
    shortest_edge::Float64
    longest_edge::Float64
end

function basis_values(geometry::CellGeometry, point_coordinates, basis_power_pairs_list; dx = 0, dy = 0)

    basis_values_vec = Basis.raw_basis(
        point_coordinates,
        geometry.centroid_coordinates,
        geometry.hx,
        geometry.hy,
        basis_power_pairs_list;
        dx = dx,
        dy = dy,
    )

    if dx == 0 && dy == 0
        basis_values_vec -= geometry.raw_means_list
        basis_values_vec[1] = 1.0
    end

    return basis_values_vec
end

function build_cell_geometry(cell_vertices_coordinate_mat, options)
    reference_quadrature_coordinates_list, quadrature_weights_1d_list = Quadrature.gauss_legendre(options.quadrature_order)

    quadrature_coordinates_mat, quadrature_weights = Quadrature.quadrilateral_quadrature(
            cell_vertices_coordinate_mat,
            reference_quadrature_coordinates_list,
            quadrature_weights_1d_list,
        )

    cell_area = sum(quadrature_weights)
    centroid_coordinates = quadrature_coordinates_mat * quadrature_weights / cell_area

    xmin_vertex_index = argmin(cell_vertices_coordinate_mat[1, :])
    xmax_vertex_index = argmax(cell_vertices_coordinate_mat[1, :])
    ymin_vertex_index = argmin(cell_vertices_coordinate_mat[2, :])
    ymax_vertex_index = argmax(cell_vertices_coordinate_mat[2, :])

    hx = cell_vertices_coordinate_mat[1, xmax_vertex_index] - cell_vertices_coordinate_mat[1, xmin_vertex_index]
    hy = cell_vertices_coordinate_mat[2, ymax_vertex_index] - cell_vertices_coordinate_mat[2, ymin_vertex_index]

    Options.check_positive(abs(hx), "zero x scale in the moving basis")
    Options.check_positive(abs(hy), "zero y scale in the moving basis")

    basis_power_pairs_list = Basis.basis_powers(options.degree)
    number_of_basis_functions = length(basis_power_pairs_list)
    number_of_quadrature_points = length(quadrature_weights)

    raw_basis_mat = zeros(number_of_basis_functions, number_of_quadrature_points)
    basis_dx_mat = zeros(number_of_basis_functions, number_of_quadrature_points)
    basis_dy_mat = zeros(number_of_basis_functions, number_of_quadrature_points)

    for quadrature_point_index in 1:number_of_quadrature_points
        quadrature_point_coordinates = quadrature_coordinates_mat[:, quadrature_point_index]

        raw_basis_mat[:, quadrature_point_index] = Basis.raw_basis(
                quadrature_point_coordinates,
                centroid_coordinates,
                hx,
                hy,
                basis_power_pairs_list,
            )

        basis_dx_mat[:, quadrature_point_index] = Basis.raw_basis(
                quadrature_point_coordinates,
                centroid_coordinates,
                hx,
                hy,
                basis_power_pairs_list;
                dx = 1,
            )

        basis_dy_mat[:, quadrature_point_index] = Basis.raw_basis(
                quadrature_point_coordinates,
                centroid_coordinates,
                hx,
                hy,
                basis_power_pairs_list;
                dy = 1,
            )
    end

    raw_means = raw_basis_mat * quadrature_weights / cell_area

    basis_mat = zeros(number_of_basis_functions, number_of_quadrature_points)

    for quadrature_point_index in 1:number_of_quadrature_points
        basis_mat[:, quadrature_point_index] = raw_basis_mat[:, quadrature_point_index] - raw_means
        basis_mat[1, quadrature_point_index] = 1.0
    end

    mass_mat = zeros(number_of_basis_functions, number_of_basis_functions)
    for quadrature_point_index in 1:number_of_quadrature_points
        basis_values_vec = basis_mat[:, quadrature_point_index]
        mass_mat += quadrature_weights[quadrature_point_index] *
            basis_values_vec * transpose(basis_values_vec)
    end

    mass_mat[1, 1] = cell_area
    mass_mat[1, 2:end] .= 0.0
    mass_mat[2:end, 1] .= 0.0

    edge_lengths = Float64[]
    for local_edge_index in 1:4
        next_local_vertex_index = mod1(local_edge_index + 1, 4)
        edge_length = norm(
            cell_vertices_coordinate_mat[:, next_local_vertex_index] - cell_vertices_coordinate_mat[:, local_edge_index],
        )
        if edge_length > 1.0e-14
            push!(edge_lengths, edge_length)
        end
    end

    scale_vertex_indices_list = [
        xmin_vertex_index,
        xmax_vertex_index,
        ymin_vertex_index,
        ymax_vertex_index,
    ]

    return CellGeometry(
        cell_vertices_coordinate_mat,
        quadrature_coordinates_mat,
        quadrature_weights,
        cell_area,
        centroid_coordinates,
        hx,
        hy,
        scale_vertex_indices_list,
        raw_means,
        mass_mat,
        basis_mat,
        basis_dx_mat,
        basis_dy_mat,
        minimum(edge_lengths),
        maximum(edge_lengths),
    )
end

function build_cells_geometries(vertices_coordinates_mat, mesh::Meshes.Mesh, options::Options.DGOptions)
    cells_geometries = [
        build_cell_geometry(
            vertices_coordinates_mat[:, mesh.cells_vertex_indices_mat[:, cell_index]],
            options,
        )
        for cell_index in 1:size(mesh.cells_vertex_indices_mat, 2)
    ]
    return cells_geometries
end

function compute_cells_coefficients(cells_moments_array, cells_geometries)
    cells_coefficients = [
        cells_geometries[cell_index].mass_mat \ cells_moments_array[:, :, cell_index]
        for cell_index in eachindex(cells_geometries)
    ]
    return cells_coefficients
end

function evaluate_solution(coefficients_mat, geometry, point_coordinates, basis_power_pairs_list; dx = 0, dy = 0)
    basis_values_vec = basis_values(geometry, point_coordinates, basis_power_pairs_list; dx = dx, dy = dy)
    state_vec = transpose(coefficients_mat) * basis_values_vec
    return state_vec
end

function positivity_sampling_points(geometry, options)
    sampling_coordinates_list = [
        geometry.quadrature_coordinates_mat[:, quadrature_point_index]
        for quadrature_point_index in eachindex(geometry.quadrature_weights_list)
    ]

    for local_vertex_index in 1:4
        push!(
            sampling_coordinates_list,
            geometry.vertex_coordinates_mat[:, local_vertex_index],
        )
    end

    edge_reference_quadrature_coordinates_list, _ = Quadrature.gauss_legendre(options.quadrature_order)
    for local_edge_index in 1:4
        next_local_vertex_index = mod1(local_edge_index + 1, 4)
        for reference_coordinate in edge_reference_quadrature_coordinates_list
            edge_coordinate = 0.5 * (1.0 + reference_coordinate)
            point_coordinates = (1.0 - edge_coordinate) *
                geometry.vertex_coordinates_mat[:, local_edge_index] + edge_coordinate *
                geometry.vertex_coordinates_mat[:, next_local_vertex_index]
            push!(sampling_coordinates_list, point_coordinates)
        end
    end

    return sampling_coordinates_list
end

end # module Geometry
