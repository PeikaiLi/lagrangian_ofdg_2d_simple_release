module Visualization

    using Plots

    import ..Geometry
    import ..Euler
    import ..States
    import ..Meshes
    import ..Discretization

    function plot_mesh(
        state::States.DGState,
        mesh::Meshes.Mesh;
        limits = nothing,
        title = "Mesh",
    )
        figure = plot(
            aspect_ratio = :equal,
            xlabel = "x",
            ylabel = "y",
            legend = false,
            title = title,
        )

        for cell_index in 1:Meshes.number_of_cells(mesh)
            coords = States.cell_vertices_coords(state, mesh, cell_index)

            plot!(
                figure,
                coords[1, [1, 2, 3, 4, 1]],
                coords[2, [1, 2, 3, 4, 1]];
                color = :black,
                linewidth = 0.5,
            )
        end

        if limits !== nothing
            plot!(figure; xlims = limits[1], ylims = limits[2])
        end

        return figure
    end

    
    function field_value(cons_vec, gamma, field::Symbol)
        if field == :density
            return cons_vec[1]
        elseif field == :pressure
            return Euler.pressure(cons_vec, gamma)
        elseif field == :velocity_x
            return cons_vec[2] / cons_vec[1]
        elseif field == :velocity_y
            return cons_vec[3] / cons_vec[1]
        elseif field == :energy
            return cons_vec[4]
        end

        throw(ArgumentError("Unknown field: $field"))
    end


    # One quadrilateral and one value per computational cell.For simplicity.
    # we can eval any pts in fact.
    function sample_field(
        state::States.DGState,
        case::Discretization.DGCase;
        field::Symbol = :density,
    )
        cells_geometries = Geometry.build_cells_geometries(
            state.vertices_coords_mat,
            case.mesh,
            case.options,
        )

        cells_coords_ls = Matrix{Float64}[]
        values_ls = Float64[]

        for cell_index in eachindex(cells_geometries)
            geometry = cells_geometries[cell_index]

            # The first basis function is 1, so its moment is the cell integral.
            mean_cons_vec =
                state.cells_moments_array[1, :, cell_index] / geometry.area

            value = field_value(
                mean_cons_vec,
                case.cells_gamma_ls[cell_index],
                field,
            )

            push!(cells_coords_ls, geometry.cell_vertices_coords_mat)
            push!(values_ls, value)
        end

        return (
            cells_coords_ls = cells_coords_ls,
            values_ls = values_ls,
        )
    end


    function color_limits(values_ls)
        lower, upper = extrema(values_ls)

        # A constant field still needs a nonzero color range.
        # if lower == upper ## plot weird!
        if isapprox(lower, upper; rtol = 1.0e-8, atol = 1.0e-12)
            padding = max(abs(lower) * 1.0e-6, 1.0e-12)
            lower -= padding
            upper += padding
        end

        return lower, upper
    end


    function draw_field(
        data;
        clims = color_limits(data.values_ls),
        limits = nothing,
        title = "Density",
        colorbar_title = "density",
    )
        figure = plot(
            aspect_ratio = :equal,
            xlabel = "x",
            ylabel = "y",
            legend = false,
            title = title,
        )

        for (coords, value) in zip(data.cells_coords_ls, data.values_ls)
            # Fill the entire cell with the color corresponding to its value.
            plot!(
                figure,
                Shape(coords[1, :], coords[2, :]);
                fill_z = [value],
                fillcolor = :viridis, # a color bar I like
                clims = clims,
                linewidth = 0,
                colorbar = true,
                colorbar_title = colorbar_title,
                label = false,
            )
        end

        if limits !== nothing
            plot!(figure; xlims = limits[1], ylims = limits[2])
        end

        return figure
    end


    function plot_field(
        state::States.DGState,
        case::Discretization.DGCase;
        field::Symbol = :density,
        time::Real = 0.0,
        clims = nothing,
        limits = nothing,
    )
        data = sample_field(state, case; field = field)

        if clims === nothing
            clims = color_limits(data.values_ls)
        end

        return draw_field(
            data;
            clims = clims,
            limits = limits,
            title = "$(field), t=$(round(time; digits=4))",
            colorbar_title = string(field),
        )
    end


    # Find a coordinate range containing the mesh at every saved time.
    function animation_limits(states_ls)
        xmin = minimum(
            minimum(state.vertices_coords_mat[1, :]) for state in states_ls
        )
        xmax = maximum(
            maximum(state.vertices_coords_mat[1, :]) for state in states_ls
        )
        ymin = minimum(
            minimum(state.vertices_coords_mat[2, :]) for state in states_ls
        )
        ymax = maximum(
            maximum(state.vertices_coords_mat[2, :]) for state in states_ls
        )

        padding = 0.02 * max(xmax - xmin, ymax - ymin, 1.0e-12)

        return (
            (xmin - padding, xmax + padding),
            (ymin - padding, ymax + padding),
        )
    end


    function draw_gif(
        solution,
        filename = "solution.gif";
        field::Symbol = :density,
        fps::Int = 10,
    )
        states_ls = solution.saved_states_ls
        times_ls = solution.saved_times_ls

        length(states_ls) == length(times_ls) ||
            throw(DimensionMismatch("States and times must match"))

        # Prepare the cell values once per saved state.
        frames_ls = [
            sample_field(state, solution.case; field = field)
            for state in states_ls
        ]

        # Keep the axes and color scale fixed throughout the animation.
        limits = animation_limits(states_ls)

        lower = minimum(minimum(data.values_ls) for data in frames_ls)
        upper = maximum(maximum(data.values_ls) for data in frames_ls)
        clims = color_limits((lower, upper))

        # refer the link below for the animation code(at bottom of the page)
        # https://persson.berkeley.edu/Programming_for_Mathematical_Applications/content/Differential_Equations/Higher_Order_Derivatives_Systems_Of_ODEs.html
        # Create an animation object
        anim = Animation()

        # Loop through every computed time step
        for frame_index in eachindex(states_ls)
            figure = draw_field(
                frames_ls[frame_index];
                clims = clims,
                limits = limits,
                title = "$(field), t=$( round(times_ls[frame_index]; digits=4) )",
                colorbar_title = string(field),
            )

            frame(anim, figure)
        end

        mkpath(dirname(abspath(filename)))

        # Compile all the frames into a GIF, specifying frames per second (fps)
        return gif(anim, filename; fps = fps)
    end

    
end # module Visualization