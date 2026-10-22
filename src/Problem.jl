module Problems

    abstract type EulerProblem end
    # Required attribute for subclass
    #   name::Symbol, smooth::Bool,....

    # methods:
    function initial_cons(problem::EulerProblem, point_coordinates)
        error("No initial state is defined for $(problem.name)")
    end

    function initial_gamma(problem::EulerProblem, point_coordinates)
        error("No initial material gamma is defined for $(problem.name)")
    end

    function exact_cons(problem::EulerProblem, point_coordinates, time)
        error("No smooth exact solution is defined for $(problem.name)")
    end

    function source_term(problem::EulerProblem, point_coordinates, time)
        return zeros(4)
    end

    has_source(problem::EulerProblem) = false

end # module Problems