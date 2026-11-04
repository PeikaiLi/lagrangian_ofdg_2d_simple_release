module Plotting

    using Plots
    import ..Basis
    import ..Geometry
    import ..Quadrature
    import ..Euler
    import ..Problems
    import ..States
    import ..Meshes
    import ..Discretization

    function plot_mesh(state::States.DGState, mesh::Meshes.Mesh;
                    limits=nothing, title="Current mesh")
        figure = plot(; aspect_ratio=:equal, xlabel="x", ylabel="y", legend=false, title=title)
        for cell_index in 1:Meshes.number_of_cells(mesh)
            coords = States.cell_vertices_coords(state, mesh, cell_index)
            plot!(figure, coords[1,[1,2,3,4,1]], coords[2,[1,2,3,4,1]];
                color=:black, linewidth=0.5)
        end
        if limits !== nothing
            plot!(figure; xlims=limits[1], ylims=limits[2])
        end
        return figure
    end
end # module Visualization
