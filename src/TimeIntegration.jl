module TimeIntegration  

    import ..DGOperator
    import ..Euler
    import ..Geometry
    import ..Problems
    import ..Discretization
    import ..Meshes
    import ..MeshMotion
    import ..Positivity
    import ..Options
    import ..States

    function stage_rhs(
        state::States.DGState,
        case::Discretization.DGCase,
        time::Real,
        )
        if Problems.has_source(case.problem)
            throw(ArgumentError("This DGOperator implements no volume source; $(case.problem.name) requires one."))
        end

        cells_geometries = Geometry.build_cells_geometries(
            state.vertices_coords_mat, case.mesh, case.options,
        )
        cells_coefficients = Geometry.moments2coefficients(
            state.cells_moments_array, cells_geometries,
        )
        vertices_velocities = MeshMotion.compute_vertex_velocities(
            state, case, cells_geometries, cells_coefficients,
        )
        cells_rhs_array = DGOperator.assemble_volume_rhs(
            state, case, time, cells_geometries, cells_coefficients, vertices_velocities,
        )
        DGOperator.assemble_flux_rhs!(
            cells_rhs_array, state, case, cells_geometries, cells_coefficients,
            vertices_velocities,
        )

        boundary_rate_ls = vec(sum(cells_rhs_array[1, :, :], dims=2)) # conservation law's test
        return (
            cells_rhs_array=cells_rhs_array,
            vertices_velocities=vertices_velocities, 
            cells_geometries = cells_geometries,
            boundary_rate_ls =boundary_rate_ls
        )
    end


    function forward_euler_step(
            state::States.DGState,
            case::Discretization.DGCase,
            time::Real,
            dt::Real,
        )
        if case.options.use_ofdg
            throw(ArgumentError("still working on OF"))
        end
        if case.options.use_positivity_limiter
            throw(ArgumentError("still working on positivity_limiter"))
        end
        if !isfinite(dt) || dt <= 0.0
            throw(ArgumentError("dt must be finite and positive"))
        end
        F_stage = stage_rhs(state, case, time)
        # cells_rhs_array, vertices_velocities, boundary_rate

        next_state = Positivity.check_stage(
            States.DGState(
                state.vertices_coords_mat + dt * F_stage.vertices_velocities,
                state.cells_moments_array + dt * F_stage.cells_rhs_array,
            ), case,
        )
        return next_state, dt * F_stage.boundary_rate_ls
    end

    
    function rk4_step(
        state::States.DGState,
        case::Discretization.DGCase,
        time::Real,
        dt::Real,
        )
        if case.options.use_ofdg
            throw(ArgumentError("still working on OF"))
        end
        if case.options.use_positivity_limiter
            throw(ArgumentError("still working on positivity_limiter"))
        end
        if !isfinite(dt) || dt <= 0.0
            throw(ArgumentError("dt must be finite and positive"))
        end
        initial_coords = state.vertices_coords_mat
        initial_moments = state.cells_moments_array

        k1_data = stage_rhs(state, case, time)

        stage2_state = Positivity.check_stage(
            States.DGState(
                initial_coords + 0.5 * dt * k1_data.vertices_velocities,
                initial_moments + 0.5 * dt * k1_data.cells_rhs_array,
            ), case,
        )
        k2_data = stage_rhs(stage2_state, case, time + 0.5 * dt)

        stage3_state = Positivity.check_stage(
            States.DGState(
                initial_coords + 0.5 * dt * k2_data.vertices_velocities,
                initial_moments + 0.5 * dt * k2_data.cells_rhs_array,
            ), case,
        )
        k3_data = stage_rhs(stage3_state, case, time + 0.5 * dt)

        stage4_state = Positivity.check_stage(
            States.DGState(
                initial_coords + dt * k3_data.vertices_velocities,
                initial_moments + dt * k3_data.cells_rhs_array,
            ), case,
        )
        k4_data = stage_rhs(stage4_state, case, time + dt)

        next_coords = initial_coords + (dt / 6.0) *
            (k1_data.vertices_velocities + 
                2.0 * k2_data.vertices_velocities + 
                2.0 * k3_data.vertices_velocities +
                k4_data.vertices_velocities)
        next_moments = initial_moments + (dt / 6.0) *
            (k1_data.cells_rhs_array + 
                2.0 * k2_data.cells_rhs_array + 
                2.0 * k3_data.cells_rhs_array +
                k4_data.cells_rhs_array)

        next_state = Positivity.check_stage(
            States.DGState(next_coords, next_moments), case,
        )

        boundary_change_ls = (dt / 6.0) *
            (k1_data.boundary_rate_ls + 
            2.0 * k2_data.boundary_rate_ls + 
            2.0 * k3_data.boundary_rate_ls + 
            k4_data.boundary_rate_ls)
        
        return next_state, boundary_change_ls
    end


    function ssprk3_step(
            state::States.DGState,
            case::Discretization.DGCase,
            time::Real,
            dt::Real,
        )

        if case.options.use_ofdg
            throw(ArgumentError("still working on OF"))
        end
        if case.options.use_positivity_limiter
            throw(ArgumentError("still working on positivity_limiter"))
        end
        if !isfinite(dt) || dt <= 0.0
            throw(ArgumentError("dt must be finite and positive"))
        end

        F_stage = stage_rhs(
            state,
            case,
            time,
            )
        
        state1 = Positivity.check_stage(
            States.DGState(
                state.vertices_coords_mat + dt * F_stage.vertices_velocities,
                state.cells_moments_array + dt * F_stage.cells_rhs_array,
            ), case,
        )

        F_stage1 = stage_rhs(
            state1,
            case,
            time + dt,
            )

        state2 = Positivity.check_stage(
            States.DGState(
                (3/4) * state.vertices_coords_mat + 
                    (1/4) * (state1.vertices_coords_mat +
                        dt * F_stage1.vertices_velocities),
                (3/4) * state.cells_moments_array + 
                    (1/4) * (state1.cells_moments_array +
                        dt * F_stage1.cells_rhs_array),
            ), case,
        )

        F_stage2 = stage_rhs(
            state2,
            case,
            time + 0.5 * dt,
            )

        state3 = Positivity.check_stage(
            States.DGState(
                (1/3) * state.vertices_coords_mat + 
                    (2/3) * (state2.vertices_coords_mat +
                        dt * F_stage2.vertices_velocities),
                (1/3) * state.cells_moments_array + 
                    (2/3) * (state2.cells_moments_array +
                        dt * F_stage2.cells_rhs_array),
            ), case,
        )
        
        boundary_change_ls =  dt * (
            1/6 * F_stage.boundary_rate_ls
            + 1/6 * F_stage1.boundary_rate_ls
            + 2/3 * F_stage2.boundary_rate_ls
        )
        return state3, boundary_change_ls
    end 


    function select_step_function(time_method::Symbol)
        if time_method == :ForwardEuler
            return forward_euler_step
        elseif time_method == :RK4
            return rk4_step
        elseif time_method == :SSPRK3
            return ssprk3_step
        end
        throw(ArgumentError("Wrong method name"))
    end


    function cfl_time_step(
        state::States.DGState,
        case::Discretization.DGCase,
        )

        mesh = case.mesh
        options = case.options
        cells_gamma_ls = case.cells_gamma_ls

        cells_geometries = Geometry.build_cells_geometries(state.vertices_coords_mat, mesh, options)
        cells_coefficients = Geometry.moments2coefficients(state.cells_moments_array, cells_geometries)

        courant_number = options.cfl

        dt = Inf
        for cell_index in eachindex(cells_geometries)
            rho, velocity_x, velocity_y, pressure_value = Euler.cons2prim(
                cells_coefficients[1, :, cell_index],
                cells_gamma_ls[cell_index],
            )

            sound_speed = sqrt(cells_gamma_ls[cell_index] * pressure_value / rho)
            length_scale = cells_geometries[cell_index].shortest_edge

            candidate_dt = courant_number * length_scale / (sound_speed + sqrt(velocity_x^2 + velocity_y^2))
            dt = min(dt, candidate_dt)
        end

    return dt
end



end # module TimeIntegration  