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
        vertices_velocities_mat = MeshMotion.compute_vertex_velocities(
            state, case, cells_geometries, cells_coefficients,
        )
        cells_rhs_array = DGOperator.assemble_volume_rhs(
            state, case, time, cells_geometries, cells_coefficients, vertices_velocities_mat,
        )
        DGOperator.assemble_flux_rhs!(
            cells_rhs_array, state, case, cells_geometries, cells_coefficients,
            vertices_velocities_mat,
        )

        boundary_rate_ls = vec(sum(cells_rhs_array[1, :, :], dims=2)) # conservation law's test
        return (
            cells_rhs_array=cells_rhs_array,
            vertices_velocities_mat=vertices_velocities_mat, 
            cells_geometries = cells_geometries,
            cells_coefficients = cells_coefficients,
            boundary_rate_ls=boundary_rate_ls,
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
        stage_data = stage_rhs(state, case, time)
        # moments_rhs, vertices_velocities, boundary_rate

        next_state = Positivity.check_stage(
            States.DGState(
                state.vertices_coords_mat + dt * stage_data.vertices_velocities,
                state.cells_moments_array + dt * stage_data.moments_rhs,
            ), case,
        )
        return next_state, dt * boundary_rate_ls
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

        k1_moments, k1_coords, k1_boundary = stage_rhs(state, case, time)

        stage2_state = Positivity.check_stage(
            States.DGState(
                initial_coords + 0.5 * dt * k1_coords,
                initial_moments + 0.5 * dt * k1_moments,
            ), case,
        )
        k2_moments, k2_coords, k2_boundary = stage_rhs(stage2_state, case, time + 0.5 * dt)

        stage3_state = Positivity.check_stage(
            States.DGState(
                initial_coords + 0.5 * dt * k2_coords,
                initial_moments + 0.5 * dt * k2_moments,
            ), case,
        )
        k3_moments, k3_coords, k3_boundary = stage_rhs(stage3_state, case, time + 0.5 * dt)

        stage4_state = Positivity.check_stage(
            States.DGState(
                initial_coords + dt * k3_coords,
                initial_moments + dt * k3_moments,
            ), case,
        )
        k4_moments, k4_coords, k4_boundary = stage_rhs(stage4_state, case, time + dt)

        next_coords = initial_coords + (dt / 6.0) *
            (k1_coords + 2.0 * k2_coords + 2.0 * k3_coords + k4_coords)
        next_moments = initial_moments + (dt / 6.0) *
            (k1_moments + 2.0 * k2_moments + 2.0 * k3_moments + k4_moments)
        next_state = check_stage(
            States.DGState(next_coords, next_moments), case,
        )

        boundary_change_ls = (dt / 6.0) *
            (k1_boundary + 2.0 * k2_boundary + 2.0 * k3_boundary + k4_boundary)
        
        return next_state, boundary_change_ls
    end


    function ssprk3_step(
            state::States.DGState,
            case::Discretization.DGCase,
            time::Real,
            dt::Real,
        )
        # finish it tmr
        # return
    end

    function select_step_function(time_method::Symbol)
        if time_method == :ForwardEuler
            return forward_euler_step
        elseif time_method == :RK4
            return rk4_step
        # elseif time_method == :SSPRK3
        #     return ssprk3_step
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