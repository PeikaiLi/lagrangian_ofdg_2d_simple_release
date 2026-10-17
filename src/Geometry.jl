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
end

function build_CellGeometry(vertex_coordinates_mat, options)
    quadrature_coordinates_list, quadrature_weights = Quadrature.quadrilateral_quadrature(
        vertex_coordinates_mat,
        options.quadrature_order,
    )


        return CellGeometry(
            vertex_coordinates_mat,
            quadrature_coordinates_list,
            quadrature_weights,
        )
end

end # module Geometry
