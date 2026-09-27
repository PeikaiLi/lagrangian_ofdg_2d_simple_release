module Types


using LinearAlgebra


struct DGOptions
    degree::Int
    use_ofdg::Bool
    use_characteristic_decomposition::Bool
    use_positivity_limiter::Bool
    cfl::Float64
    quadrature_order::Int
    density_floor::Float64
    pressure_floor::Float64
    max_retries::Int
    max_steps::Int
end


struct Mesh
    initial_boundary_box::Vector{Float64} # x0,x1,y0,y1
    nx::Int # number of cells at x-coordinate
    ny::Int
    cells_vertex_indices_mat::Matrix{Int}
    faces_list::Vector{Face}
    
    periodic::Bool
    vertex2group_index::Vector{Int}
    
    polar::Bool
    has_origin::Bool
end

struct DGState
    # 2 × number_of_vertices. Column `vertex_index` is one physical vertex.
    vertices_coordinates_mat::Matrix{Float64}
end

struct Face

    start_vertex_index::Int
    end_vertex_index::Int

    # respect to the direction of the face,
    # left cell is on the left side of the face, right cell is on the right side of the face
    left_cell_index::Int 
    right_cell_index::Int

    contact::Bool

    boundary::Symbol
    periodic_shift_vec::Vector{Float64}    
end