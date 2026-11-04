module Diagnostics

    import ..States
    import ..Geometry
    import ..Quadrature
    import ..Basis
    import ..Problems
    import ..Discretization
    import ..Solver

    function total_cons(state::States.DGState)
        return vec(sum(state.cells_moments_array[1, :, :], dims=2))
    end

    function conservation_report(solution)
        initial_total_ls = total_cons(solution.case.initial)
        final_total_ls = total_cons(solution.state)
        change_ls = final_total_ls - initial_total_ls
        return (
            initial_total_ls = initial_total_ls,
            final_total_ls = final_total_ls,
            change_ls = change_ls,
            boundary_change_ls = solution.boundary_change_ls,
            balance_error_ls = change_ls - solution.boundary_change_ls, # It must be zero.
        )
    end

end # module Diagnostics