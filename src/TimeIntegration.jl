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

end # module TimeIntegration  