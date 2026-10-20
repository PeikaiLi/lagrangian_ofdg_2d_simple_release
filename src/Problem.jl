module Problems

    abstract type EulerProblem end
    # Required attribute: name::Symbol, smooth::Bool, final_time::Float64

    # methods:
    function initial_conserve_vec(problem::EulerProblem, point_coordinates)
        error("No initial state is defined for $(problem.name)")
    end

    function initial_gamma(problem::EulerProblem, point_coordinates)
        error("No initial material gamma is defined for $(problem.name)")
    end

    function exact_state(problem::EulerProblem, point_coordinates, time)
        error("No smooth exact solution is defined for $(problem.name)")
    end

    function source_term(problem::EulerProblem, point_coordinates, time)
        return zeros(4)
    end

    has_source(problem::EulerProblem) = false
    radial_profile(problem::EulerProblem) = false

end # module Problems