module Positivity

import ..Basis
import ..Euler
import ..Geometry
import ..States
import ..Meshes
import ..Options

function check_state(
    state::States.DGState,
    case::Discretization.DGCase
    )
    mesh = case.mesh
    cells_gamma_ls = case.cells_gamma_ls
    options = case.options
    cells_geometries = Geometry.build_cells_geometries(state.vertices_coords_mat, mesh, options)
    cells_coefficients = Geometry.moments2coefficients(state.cells_moments_array, cells_geometries)
    exponents_ls = Basis.basis_powers(options.degree)

    for cell_index in eachindex(cells_geometries)
        cell_average_cons_vec = cells_coefficients[1, :, cell_index]
        gamma = cells_gamma_ls[cell_index]

        if !Euler.is_admissible(cell_average_cons_vec, gamma, options)
            throw(DomainError(cell_index,
                "nonphysical cell average: rho=$(cell_average_cons_vec[1]), " *
                "p=$(Euler.pressure(cell_average_cons_vec, gamma))"))
        end

        for point_coord in cells_geometries[cell_index].sampling_coords_ls
            cons_vec = Geometry.eval_cons(
                cells_coefficients[:, :, cell_index],
                cells_geometries[cell_index],
                point_coord,
                exponents_ls,
            )
            if !Euler.is_admissible(cons_vec, gamma, options)
                throw(DomainError(cell_index,
                    "nonphysical polynomial value at point=$(Tuple(point_coord)): " *
                    "rho=$(cons_vec[1]), p=$(Euler.pressure(cons_vec, gamma))"))
            end
        end
    end

    return state
end

end # module Positivity
