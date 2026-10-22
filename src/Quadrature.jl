module Quadrature

    using LinearAlgebra: diagm, eigen

    function orient2d_det(a, b, c)
        # >0 if triabc is counter-clockwise, <0 if clockwise, 0 if collinear
        # a, b, c are 2D points (vectors of length 2)
        return ((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
    end


    # Below, refer to UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>
    # https://nbviewer.org/url/persson.berkeley.edu/math228b/notebooks-plots/high_order_utils.ipynb
    """
        function gauss_quad(p)

    Gaussian quadrature on [-1,1] for given degree of precision `p`
    """
    function gauss_quad(p)
        n = ceil(Int,(p+1)/2)
        b = 1:n-1
        b = @. b / sqrt(4*b^2 - 1)
        eval, evec = eigen(diagm(1 => b, -1 => b))
        return eval, 2*evec[1,:].^2
    end


    function quadrilateral_quadrature(cell_vertices_coord_mat,quadrature_order)
    # check mesh is valid quadrilateral
        or_2area_triangle123 = orient2d_det(
            cell_vertices_coord_mat[:, 1],
            cell_vertices_coord_mat[:, 2],
            cell_vertices_coord_mat[:, 3],
        )
        or_2area_triangle134 = orient2d_det(
            cell_vertices_coord_mat[:, 1],
            cell_vertices_coord_mat[:, 3],
            cell_vertices_coord_mat[:, 4],
        )
        or_2area_triangle124 = orient2d_det(
            cell_vertices_coord_mat[:, 1],
            cell_vertices_coord_mat[:, 2],
            cell_vertices_coord_mat[:, 4],
        )
        or_2area_triangle234 = orient2d_det(
            cell_vertices_coord_mat[:, 2],
            cell_vertices_coord_mat[:, 3],
            cell_vertices_coord_mat[:, 4],
        )


        if or_2area_triangle123 >= 0 && or_2area_triangle134 >=0
            two_triangle_vertex_indices_list = [(1, 2, 3), (1, 3, 4)]
        elseif or_2area_triangle124 < 0 || or_2area_triangle234 < 0
            throw(DomainError("crossed or inverted quadrilateral"))
        else
            two_triangle_vertex_indices_list = [(1, 2, 4), (2, 3, 4)]
        end

        # quadrature points and weights for this quadrilateral
        quadrature1d_pts_ls, quadrature1d_wts_ls = gauss_quad(quadrature_order)

        # scale to [-1,1] -> [0,1]
        quadrature1d_pts_ls = quadrature1d_pts_ls * 0.5 .+ 0.5 
        quadrature1d_wts_ls = quadrature1d_wts_ls * 0.5 
    
        quadrilateral_quadrature_pts_ls = Vector{Float64}[]
        quadrilateral_quadrature_wts_ls = Float64[]
        for (v1, v2, v3) in two_triangle_vertex_indices_list
            a = cell_vertices_coord_mat[:, v1]
            b = cell_vertices_coord_mat[:, v2]
            c = cell_vertices_coord_mat[:, v3]
            X(alpha, beta) = a + alpha * (b - a) + (1-alpha)*beta * (c - a)
            for i in 1:length(quadrature1d_pts_ls)
                for j in 1:length(quadrature1d_pts_ls)
                    alpha = quadrature1d_pts_ls[i]
                    beta = quadrature1d_pts_ls[j]
                    weight = (1-alpha)*quadrature1d_wts_ls[i] * quadrature1d_wts_ls[j] * abs(orient2d_det(a, b, c))

                    push!(quadrilateral_quadrature_pts_ls, X(alpha, beta))
                    push!(quadrilateral_quadrature_wts_ls, weight)
                end
            end
        end

        cell_area = sum(quadrilateral_quadrature_wts_ls)
        if !isfinite(cell_area) || cell_area <= 0.0
            throw(DomainError(cell_area, "invalid cell area"))
        end

        return quadrilateral_quadrature_pts_ls, quadrilateral_quadrature_wts_ls
    end
end
