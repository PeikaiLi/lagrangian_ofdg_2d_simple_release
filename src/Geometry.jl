module Geometry

using LinearAlgebra: norm
import ..Basis
import ..Meshes
import ..Options
import ..Quadrature

struct CellGeometry
    vertex_coordinates_mat::Matrix{Float64}
    quadrature_coordinates_ls::Vector{Vector{Float64}} 
    quadrature_weights_ls::Vector{Float64}
    area::Float64
    centroid_coordinates::Vector{Float64}
    hx::Float64
    hy::Float64
end

function build_CellGeometry(vertex_coordinates_mat, options)

    quadrature_coordinates_list, quadrature_weights = Quadrature.quadrilateral_quadrature(
        vertex_coordinates_mat,
        options.quadrature_order,
    )

    area = sum(quadrature_weights)

    centroid_coordinates = zeros(2)
    for i in length(quadrature_coordinates_list)
        centroid_coordinates += quadrature_weights[i] *
            quadrature_coordinates_list[i]
    end
    centroid_coordinates = centroid_coordinates/area
       
    xmin_vertex_index = argmin(vertex_coordinates_mat[1, :])
    xmax_vertex_index = argmax(vertex_coordinates_mat[1, :])
    ymin_vertex_index = argmin(vertex_coordinates_mat[2, :])
    ymax_vertex_index = argmax(vertex_coordinates_mat[2, :])

    hx = vertex_coordinates_mat[1, xmax_vertex_index] - vertex_coordinates_mat[1, xmin_vertex_index]
    hy = vertex_coordinates_mat[2, ymax_vertex_index] - vertex_coordinates_mat[2, ymin_vertex_index]





        return CellGeometry(
            vertex_coordinates_mat,
            quadrature_coordinates_list,
            quadrature_weights,
            area,
            centroid_coordinates,
            hx,
            hy


        )
end

end # module Geometry
