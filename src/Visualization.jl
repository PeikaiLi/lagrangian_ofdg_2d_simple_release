module Visualization

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


    function field_value(cons_vec, gamma, field::Symbol)
        field == :density && return cons_vec[1]
        field == :energy && return cons_vec[4]
        field == :pressure && return Euler.pressure(cons_vec, gamma)
        field == :velocity_x && return cons_vec[2]/cons_vec[1]
        field == :velocity_y && return cons_vec[3]/cons_vec[1]
        throw(ArgumentError("wrong name"))
    end


    function sample_field(state::States.DGState, case::Discretization.DGCase;
                        subdivisions::Int=3, field::Symbol=:density,
                        mode::Symbol=:numerical, time::Real=0.0)
    
        mode in (:numerical,:exact,:error) || throw(ArgumentError("unknown display mode"))

        geometries = Geometry.build_cells_geometries(state.vertices_coords_mat,case.mesh,case.options)
        coefficients = Geometry.moments2coefficients(state.cells_moments_array,geometries)
        exponents = Basis.basis_powers(case.options.degree)

        triangles_ls = Matrix{Float64}[]
        values_ls = Float64[]
        for cell_index in eachindex(geometries)
            geometry = geometries[cell_index]
            coords = geometry.cell_vertices_coords_mat
            # Use the same valid diagonal convention as the cell quadrature.
            if Quadrature.orient2d_det(coords[:,1],coords[:,2],coords[:,3]) >= 0 &&
            Quadrature.orient2d_det(coords[:,1],coords[:,3],coords[:,4]) >= 0
                triangle_indices = ((1,2,3),(1,3,4))
            else
                triangle_indices = ((1,2,4),(2,3,4))
            end
            for (ia,ib,ic) in triangle_indices
                a,b,c = coords[:,ia],coords[:,ib],coords[:,ic]
                if Quadrature.orient2d_det(a,b,c) == 0.0 
                     continue
                end
                point(i,j) = a + (i/subdivisions)*(b-a) + (j/subdivisions)*(c-a)
                for i in 0:subdivisions-1, j in 0:subdivisions-1-i
                    patches = [hcat(point(i,j),point(i+1,j),point(i,j+1))]
                    if i+j <= subdivisions-2
                        push!(patches,hcat(point(i+1,j),point(i+1,j+1),point(i,j+1)))
                    end
                    for triangle in patches
                        point_coord = vec(sum(triangle,dims=2))/3
                        gamma = case.cells_gamma_ls[cell_index]
                        numerical = field_value(Geometry.eval_cons(coefficients[:,:,cell_index],
                            geometry,point_coord,exponents),gamma,field)
                        value = numerical
                        if mode != :numerical
                            exact = field_value(Problems.exact_cons(case.problem,point_coord,time),gamma,field)
                            value = mode == :exact ? exact : numerical-exact
                        end
                        isfinite(value) || throw(DomainError(value,"nonfinite plotting value"))
                        push!(triangles_ls,triangle)
                        push!(values_ls,value)
                    end
                end
            end
        end
        return (triangles_ls=triangles_ls,values_ls=values_ls)
    end


    function color_limits(values_ls)
        isempty(values_ls) && throw(ArgumentError("no values to plot"))
        lo,hi = extrema(values_ls)
        if lo == hi
            padding = max(abs(lo)*1e-6,1e-12)
            lo -= padding
            hi += padding
        end
        return (lo,hi)
    end

    
    function draw_field(data; limits=nothing, palette=:viridis,
                        title="Density", colorbar_title="density")

        clims = color_limits(data.values_ls) 
        (all(isfinite,clims) && clims[1] < clims[2]) || throw(ArgumentError("invalid color limits"))

        colors = cgrad(palette)
        figure = plot(;aspect_ratio=:equal,xlabel="x",ylabel="y",legend=false,title=title)
        for (triangle,value) in zip(data.triangles_ls,data.values_ls)
            fraction = clamp((value-clims[1])/(clims[2]-clims[1]),0.0,1.0)
            plot!(figure,Shape(triangle[1,:],triangle[2,:]);
                fillcolor=colors[fraction],linecolor=colors[fraction],linewidth=0,label=false)
        end
        # Invisible markers provide a colorbar matching the explicit polygon colors.
        point_coord = data.triangles_ls[1][:,1]
        scatter!(figure,fill(point_coord[1],2),fill(point_coord[2],2);
                marker_z=collect(clims),color=palette,clims=clims,markersize=0,
                markerstrokewidth=0,colorbar=true,colorbar_title=colorbar_title,label=false)
        if limits !== nothing
            plot!(figure;xlims=limits[1],ylims=limits[2])
        end
        return figure
    end






end # module Visualization
