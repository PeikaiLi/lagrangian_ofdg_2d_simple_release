module Solver
    import ..Discretization
    import ..Positivity
    import ..Problems
    import ..States
    import ..TimeIntegration
    using LinearAlgebra: SingularException


    function raw_solve_2d(case::Discretization.DGCase;
                  time_method::Symbol = :RK4,
                  dt_max::Real = Inf,
                  save_every::Int = 0,
                  verbose::Bool = false)
    
        options = case.options
        step_function = TimeIntegration.select_step_function(time_method)
        state = States.DGState(copy(case.initial.vertices_coords_mat), copy(case.initial.cells_moments_array))
        Positivity.check_state(state, case)
        time = 0.0
        number_of_steps = 0
        number_of_retries = 0
        boundary_change_ls = zeros(4)
        saved_times_ls = [time]
        saved_states_ls = [state]

        while time < case.final_time
            number_of_steps < options.max_steps || error("max_steps reached at time=$time")
            remaining_time = case.final_time - time
            cfl_dt = TimeIntegration.cfl_time_step(state, case)
            (isfinite(cfl_dt) && cfl_dt > 0.0) || throw(DomainError(cfl_dt, "invalid CFL time step"))
            dt = min(cfl_dt, Float64(dt_max), remaining_time)
            retries = 0

            # A rejected trial does not change the accepted state, time, or budget.
            while true
                try
                    next_state, step_boundary_change_ls = step_function(state, case, time, dt)
                    state = next_state
                    boundary_change_ls += step_boundary_change_ls
                    break
                catch exception
                    # Do not turn missing methods or programming errors into retries.
                    if !(exception isa DomainError || exception isa SingularException)
                        rethrow()
                    end
                    retries < options.max_retries || rethrow()
                    dt *= 0.5
                    retries += 1
                end
            end

            time = dt == remaining_time ? case.final_time : time + dt
            number_of_steps += 1
            number_of_retries += retries
            if save_every > 0 && number_of_steps % save_every == 0
                push!(saved_times_ls, time)
                push!(saved_states_ls, state)
            end
            if verbose
                println("step=", number_of_steps, " time=", time, " dt=", dt, " retries=", retries)
            end
        end

        if saved_times_ls[end] != time
            push!(saved_times_ls, time)
            push!(saved_states_ls, state)
        end

        return (
            state = state,
            case = case,
            time = time,
            time_method = time_method,
            number_of_steps = number_of_steps,
            number_of_retries = number_of_retries,
            boundary_change_ls = boundary_change_ls,
            saved_times_ls = saved_times_ls,
            saved_states_ls = saved_states_ls,
        )
        end
end