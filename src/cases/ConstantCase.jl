module ConstantCase

import ..Euler
import ..Meshes
import ..Options
import ..Problems
import ..States

# just play with a simple one.
# I need test my shit code.
struct Problem <: Problems.EulerProblem
    name::Symbol
    cells_gamma_list::Vector{Float64}
    smooth::Bool
    final_time::Float64
    nx::Int
    ny::Int
    domain_bounds_list::Vector{Float64}
end

function Problems.exact_state(problem::Problem, point_coordinates, time)
    return Problems.initial_state(problem, point_coordinates, 1)
end

function Problems.initial_state(problem::Problem, point_coordinates, cell_index)
    gamma = problem.cells_gamma_list[cell_index]
    return Euler.prim_to_cons(1.0, 0.3, -0.2, 1.0, gamma)
end

function make_case(; nx = nothing, ny = nothing, degree = 2,
                   OF = true, characteristic = true, positivity = false, cfl = nothing)
    nx = nx === nothing ? 3 : nx
    ny = ny === nothing ? nx : ny
    box = (0.0, 1.0, 0.0, 1.0)
    boundaries = (:wall, :wall, :wall, :wall)
    mesh, initial_coordinates = Meshes.make_mesh(box, nx, ny;
        boundaries = boundaries, periodic = true, polar = false,
        saltzman = false)
    cells_gamma_list = fill(1.4, nx * ny)
    problem = Problem(:constant, cells_gamma_list, true, 0.02,
                      nx, ny, collect(box))
    options = Options.DGOptions(degree = degree, OF = OF, characteristic = characteristic,
                                positivity = positivity, cfl = cfl)
    initial_function = (point_coordinates, cell_index) ->
        Problems.initial_state(problem, point_coordinates, cell_index)
    initial = States.project_initial_condition(initial_coordinates, mesh, options, initial_function)
    return Problems.DGCase(initial, mesh, problem, options)
end

end # module ConstantCase
