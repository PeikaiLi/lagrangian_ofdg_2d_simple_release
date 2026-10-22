module ConstantCase

    import ..Euler
    import ..Meshes
    import ..Options
    import ..Problems
    import ..Discretization

    struct ConstantProblem <: Problems.EulerProblem
        name::Symbol
        gamma::Float64
        smooth::Bool
        domain_bounds_list::Vector{Float64}
        boundaries::NTuple{4,Symbol}
        density::Float64
        velocity_x::Float64
        velocity_y::Float64
        pressure_value::Float64
    end

    function ConstantProblem(; density = 1.0, velocity_x = 0.3, velocity_y = -0.2,
                    pressure_value = 1.0, gamma = 1.4,
                    box = (0.0, 1.0, 0.0, 1.0), boundaries = (:periodic, :periodic, :periodic, :periodic))
        return ConstantProblem(:constant, gamma, true, collect(box), boundaries, density, velocity_x, velocity_y, pressure_value)
    end

    function Problems.initial_cons(problem::ConstantProblem, point_coordinates)
        gamma = problem.gamma
        return Euler.prim2cons(problem.density, problem.velocity_x, problem.velocity_y, problem.pressure_value, gamma)
    end

    Problems.initial_gamma(problem::ConstantProblem, point_coordinates) = problem.gamma

    function Problems.exact_cons(problem::ConstantProblem, point_coordinates, time)
        return Problems.initial_cons(problem, point_coordinates)
    end


    function make_case(problem::ConstantProblem = ConstantProblem(); nx = 3, ny = 3, degree = 2,
                    OF = true, characteristic = true, positivity = false, cfl = nothing,
                    final_time = 0.02)
        box = problem.domain_bounds_list
        boundaries = problem.boundaries
        mesh, initial_coords = Meshes.make_mesh(box, nx, ny;
            boundaries = boundaries)
        options = Options.DGOptions(degree = degree, OF = OF, characteristic = characteristic,
                                    positivity = positivity, cfl = cfl)
        return Discretization.DGCase(problem, mesh, initial_coords, options;
                                    final_time = final_time)
    end

end # module ConstantCase
