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

        extreme_vertex_indices_ls::Vector{Int}
        hx::Float64
        hy::Float64
        

        raw_basis_means_ls::Vector{Float64}
        basis_at_quadrature_mat::Matrix{Float64}
        basis_dx_at_quadrature_mat::Matrix{Float64}
        basis_dy_at_quadrature_mat::Matrix{Float64}

        mass_mat::Matrix{Float64}
        shortest_edge::Float64
        longest_edge::Float64
    end

    function build_CellGeometry(vertex_coordinates_mat, options::Options.DGOptions)

        quadrature_coordinates_list, quadrature_weights = Quadrature.quadrilateral_quadrature(
            vertex_coordinates_mat,
            options.quadrature_order,
        )
        number_of_quadrature_pts = length(quadrature_weights)

        area = sum(quadrature_weights)
        centroid_coordinates = zeros(2)
        for i in 1:length(number_of_quadrature_pts)
            centroid_coordinates += quadrature_weights[i] *
                quadrature_coordinates_list[i]
        end
        centroid_coordinates = centroid_coordinates/area
        
        xmin_vertex_index = argmin(vertex_coordinates_mat[1, :])
        xmax_vertex_index = argmax(vertex_coordinates_mat[1, :])
        ymin_vertex_index = argmin(vertex_coordinates_mat[2, :])
        ymax_vertex_index = argmax(vertex_coordinates_mat[2, :])

        extreme_vertex_indices_ls = [
            xmin_vertex_index,
            xmax_vertex_index,
            ymin_vertex_index,
            ymax_vertex_index,
        ]
        hx = vertex_coordinates_mat[1, xmax_vertex_index] - vertex_coordinates_mat[1, xmin_vertex_index]
        hy = vertex_coordinates_mat[2, ymax_vertex_index] - vertex_coordinates_mat[2, ymin_vertex_index]
        Options.check_positive(abs(hx), "zero hx scale in the moving basis")
        Options.check_positive(abs(hy), "zero hy scale in the moving basis")

        exponents_ls = Basis.basis_powers(options.degree)
        number_of_basis  = length(exponents_ls)
        number_of_quadrature_pts = length(quadrature_weights)

        raw_basis_mat = zeros(number_of_basis, number_of_quadrature_pts)
        basis_dx_mat = zeros(number_of_basis, number_of_quadrature_pts)
        basis_dy_mat = zeros(number_of_basis, number_of_quadrature_pts)

        for i in 1:length(number_of_quadrature_pts)
            pt = quadrature_coordinates_list[i]
            raw_basis_mat[:, i] = Basis.raw_basis(
                    pt,
                    centroid_coordinates,
                    hx,
                    hy,
                    basis_power_pairs_list,
                )

            basis_dx_mat[:, i] = Basis.raw_basis(
                    pt,
                    centroid_coordinates,
                    hx,
                    hy,
                    basis_power_pairs_list;
                    dx = 1,
                )

            basis_dy_mat[:, i] = Basis.raw_basis(
                    pt,
                    centroid_coordinates,
                    hx,
                    hy,
                    basis_power_pairs_list;
                    dy = 1,
                )
        end

        raw_means = raw_basis_mat * quadrature_weights / cell_area
        basis_mat = zeros(size(raw_basis_mat))
        for i in 1:length(number_of_quadrature_pts)
            basis_mat[:, i] = raw_basis_mat[:, i] - raw_means
        end
        @assert maximum(abs(basis_mat[1, :] - 1.0)) < options.tol
        basis_mat[1, :] .= 1.0

        mass_mat = zeros(number_of_basis,number_of_basis)
        for i in 1:number_of_quadrature_pts
            mass_mat += quadrature_weights*basis_mat(:,i)*basis_mat(:,i)'
        end

        # maybe I will put a test later
        mass_mat[1, 1] = cell_area
        mass_mat[1, 2:end] .= 0.0
        mass_mat[2:end, 1] .= 0.0

        edge_lengths = Float64[]
        for edge_index in 1:4
            next_index = mod1(edge_index + 1, 4)
            edge_length = norm(
                vertex_coordinates_mat[:, next_index] - vertex_coordinates_mat[:, edge_index],
            )
            if edge_length > options.tol
                push!(edge_lengths, edge_length)
            end
        end

        return CellGeometry(
            vertex_coordinates_mat,
            quadrature_coordinates_list,
            quadrature_weights,
            area,
            centroid_coordinates,
            extreme_vertex_indices_ls,
            hx,
            hy,
            raw_means,
            basis_mat,
            basis_dx_mat,
            basis_dy_mat,
            mass_mat,
            minimum(edge_lengths),
            maximum(edge_lengths),
        )
    end
    
    # Hopefully vertices_coordinates_mat and vertex_coordinates_mat is different enough
    # I will forget it one day <_<!
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



    function eval_basis_values(geometry::CellGeometry, point_coordinates, exponents_ls; dx = 0, dy = 0)
        basis_values_vec = Basis.raw_basis(
            point_coordinates,
            geometry.centroid_coordinates,
            geometry.hx,
            geometry.hy,
            exponents_ls;
            dx = dx,
            dy = dy,
        )

        if dx == 0 && dy == 0
            basis_values_vec -= geometry.raw_basis_means_ls
            basis_values_vec[1] = 1.0
        end

        return basis_values_vec
    end


end # module Geometry
