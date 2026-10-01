module Problems

import ..Meshes
import ..Options
import ..States

# Concrete problem types and their methods live in code/cases/.
abstract type EulerProblem end


# {P} declares a type parameter, used below in problem::P.
# P <: EulerProblem restricts P to EulerProblem or its subtypes.
struct DGCase{P<:EulerProblem}
    initial::States.DGState
    mesh::Meshes.Mesh
    problem::P
    options::Options.DGOptions
end

# Required method: initial_state(problem, point_coordinates, cell_index).
# Declare the shared function so each problem type can provide its own initial state.
function initial_state end

function exact_state(problem::EulerProblem, point_coordinates, time)
    error("No smooth exact solution is defined for $(problem.name)")
end

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

function boundary_pressure(problem::EulerProblem, boundary, time)
    error("No prescribed pressure is defined for $(problem.name)")
end

function time_step_settings(problem::EulerProblem, options, time)
    return options.cfl, Inf
end

end # module Problems
