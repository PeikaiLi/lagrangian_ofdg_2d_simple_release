module Diagnostics

    using Plots
    import ..States
    import ..Geometry
    import ..Quadrature
    import ..Basis
    import ..Problems
    import ..Discretization
    import ..Solver


    function total_cons(state::States.DGState)
        return vec(sum(state.cells_moments_array[1, :, :], dims=2))
    end


    function conservation_report(solution)
        initial_total_ls = total_cons(solution.case.initial)
        final_total_ls = total_cons(solution.state)
        change_ls = final_total_ls - initial_total_ls
        return (
            initial_total_ls = initial_total_ls,
            final_total_ls = final_total_ls,
            change_ls = change_ls,
            boundary_change_ls = solution.boundary_change_ls,
            balance_error_ls = change_ls - solution.boundary_change_ls, # It must be zero.
        )
    end


    function smooth_error_norms(state::States.DGState, case::Discretization.DGCase, time::Real;
                     quadrature_order::Int = max(case.options.quadrature_order, 2*case.options.degree + 5))

        case.problem.smooth || throw(ArgumentError("This error report requires a smooth exact solution"))

        cells_geometries = Geometry.build_cells_geometries(state.vertices_coords_mat, case.mesh, case.options)
        cells_coefficients = Geometry.moments2coefficients(state.cells_moments_array, cells_geometries)
        exponents_ls = Basis.basis_powers(case.options.degree)
        absolute_error_integrals_ls = zeros(4)
        squared_error_integrals_ls = zeros(4)
        sampled_Linf_ls = zeros(4)

        total_area = 0.0
        for cell_index in eachindex(cells_geometries)
            geometry = cells_geometries[cell_index]
            points_ls, weights_ls = Quadrature.quadrilateral_quadrature(geometry.cell_vertices_coords_mat, quadrature_order)
            for point_index in eachindex(weights_ls)
                point_coord = points_ls[point_index]
                cons_vec = Geometry.eval_cons(cells_coefficients[:, :, cell_index], geometry, point_coord, exponents_ls)
                exact_cons_vec = Problems.exact_cons(case.problem, point_coord, time)
                error_vec = abs.(cons_vec - exact_cons_vec)
                absolute_error_integrals_ls += weights_ls[point_index] * error_vec
                squared_error_integrals_ls += weights_ls[point_index] * error_vec.^2
                sampled_Linf_ls = max.(sampled_Linf_ls, error_vec)
            end

            # L_inf can consider vex and boundary
            for point_coord in geometry.sampling_coords_ls # repeat a lot but I just do it...
                cons_vec = Geometry.eval_cons(cells_coefficients[:, :, cell_index], geometry, point_coord, exponents_ls)
                sampled_Linf_ls = max.(sampled_Linf_ls, abs.(cons_vec - Problems.exact_cons(case.problem, point_coord, time)))
            end
            total_area += geometry.area
        end
        return (
            L1_ls = absolute_error_integrals_ls / total_area,
            L2_ls = sqrt.(squared_error_integrals_ls / total_area),
            Linf_ls = sampled_Linf_ls,
            total_area = total_area,
        )
    end

    ### need to define case_builder later
    function convergence_study_data(case_builder; mesh_sizes_ls=[4, 8, 16], degree=2,
                            final_time=0.1, time_method=:RK4, dt_max=Inf)

        rows_ls = NamedTuple[]

        for nx in mesh_sizes_ls
            case = case_builder(nx, degree, final_time) # case::Discretization.DGCase
            solution = Solver.raw_solve_2d(case; time_method=time_method, dt_max=dt_max)
            errors = smooth_error_norms(solution.state, case, solution.time)
            order_ls = fill(NaN, 4)
            if !isempty(rows_ls)
                previous = rows_ls[end]
                for component_index in 1:4
                    if previous.L2_ls[component_index] > 0 && errors.L2_ls[component_index] > 0
                        order_ls[component_index] = log(previous.L2_ls[component_index] / errors.L2_ls[component_index]) / log(nx / previous.nx)
                    end
                end
            end

            # first one just push 
            push!(rows_ls, (nx=nx, degree=degree, L1_ls=errors.L1_ls, L2_ls=errors.L2_ls,
                        Linf_ls=errors.Linf_ls, order_ls=order_ls,
                        number_of_steps=solution.number_of_steps, number_of_retries=solution.number_of_retries))
        end

        return rows_ls
    end


    function convergence_study_plot(converge_data)
        rows_ls = converge_data
        degree = rows_ls[1].degree
        # Plot the L2 err
        inverse_n_ls = [1.0 / row.nx for row in rows_ls]
        labels_ls = ["rho", "rho*u", "rho*v", "E"]

        figure = plot(
            xscale = :log10,
            yscale = :log10,
            xlabel = "1/N",
            ylabel = "L2 error",
            title = "L2 error, degree = $degree",
            legend = :best,
        )

        for component_index in 1:4
            errors_ls = [row.L2_ls[component_index] for row in rows_ls]
            plot!(
                figure,
                inverse_n_ls,
                errors_ls;
                label = labels_ls[component_index],
                marker = :circle,
                linewidth = 2,
            )
        end


        # Anchor line
        ref_index = argmax(inverse_n_ls)
        ref_x = inverse_n_ls[ref_index]
        ref_y = 0.1 * minimum(rows_ls[ref_index].L2_ls)

        # On log-log axes, y = C * (1/N)^p has slope p.
        # change a constant to fix it
        # y = ref_y .* (1/N  ./ (ref_x) )^p  
        for ref_order in 1:4
            ref_errors_ls =
                ref_y .* (inverse_n_ls ./ ref_x) .^ ref_order

            plot!(
                figure,
                inverse_n_ls,
                ref_errors_ls;
                label = "slope $ref_order",
                linestyle = :dash,
                linewidth = 1.5,
            )
        end
        display(figure)

        return rows_ls
    end
    # convergence_study_plot(converge_data)


end # module Diagnostics