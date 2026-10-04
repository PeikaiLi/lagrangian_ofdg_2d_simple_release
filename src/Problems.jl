module Problems

import ..Meshes
import ..Options
import ..States

# Concrete problem types and their methods live in code/cases/.
abstract type EulerProblem end


struct DGCase{P<:EulerProblem}
    initial::States.DGState
    mesh::Meshes.Mesh
    problem::P
    options::Options.DGOptions
end

# Most function will be defined in the concrete problem types, 
# but I will define some generic functions here.
# So that the I can call them without knowing the concrete problem type.

# Each concrete problem must define its own initial state 
# So I just write error.
function initial_state(problem::EulerProblem, point_coordinates, cell_index)
    error("No initial state is defined for $(problem.name)")
    # return the initial state[rou, rou*u, rou*v, E].
    # When I say state i always refer to Conservation Quantities. 
    # I need remind myself everywhere If I think I feel I might confused in future.
end

# I am not intending to define a smooth exact solution for any of the problems, 
# So if you don't have real solution, just do not call this fucntion.
# If you use this function, you don't define the real exaact solution, you fucked up. SO I write error.
function exact_state(problem::EulerProblem, point_coordinates, time)
    error("No smooth exact solution is defined for $(problem.name)")
    # return Conservation quantities [rou, rou*u, rou*v, E] just like initial_state.
end

# I may define a reference_state,since Shu use high order WENO as exact sol, but I can just use the exact_state
# as a reference_state, so I will not define a reference_state for now.

function source_term(problem::EulerProblem, point_coordinates, time)
    return zeros(4)
end

has_source(problem::EulerProblem) = false
radial_profile(problem::EulerProblem) = false

function wall_velocity(problem::EulerProblem, boundary, point_coordinates, time)
    if boundary in (:wall, :circular_wall)
        return [0.0, 0.0]
    end
    error("Define wall_velocity for this problem and moving-wall boundary")
end

# I will only use it for piston problems, 
# Actually, I am not sure it is a good idea defined it here, but I will leave it for now.
function boundary_pressure(problem::EulerProblem, boundary, time)
    error("No prescribed pressure is defined for $(problem.name)")
end

function time_step_settings(problem::EulerProblem, options, time)
    return options.cfl, Inf
end

end # module Problems