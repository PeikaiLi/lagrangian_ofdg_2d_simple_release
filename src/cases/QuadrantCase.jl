module QuadrantCase

    import ..Euler
    import ..Problems
    import ..Meshes
    import ..Options
    import ..Discretization

    struct QuadrantProblem <: Problems.EulerProblem
        name::Symbol
        smooth::Bool
        gamma::Float64
    end

    QuadrantProblem() = QuadrantProblem(:quadrant, false, 1.4)

    Problems.initial_gamma(problem::QuadrantProblem, point_coord) = problem.gamma


    function Problems.initial_cons(problem::QuadrantProblem, point_coord)
        x, y = point_coord
        if x < 0.5 && y < 0.5
            rho, pressure = 1.0, 1.0
        elseif x >= 0.5 && y < 0.5
            rho, pressure = 0.7, 0.4
        elseif x < 0.5 && y >= 0.5
            rho, pressure = 0.5, 0.3
        else
            rho, pressure = 0.125, 0.1
        end
        return Euler.prim2cons(rho, 0.0, 0.0, pressure, problem.gamma)
    end

    
    function make_case(; nx=20, ny=nx, degree=2, final_time=0.05,
                    OF=false, characteristic=false, positivity=false, 
                    cfl=0.03,
                    max_retries = 5
                    )

        (iseven(nx) && iseven(ny)) || throw(ArgumentError("Use even nx and ny to align the interfaces"))

        problem = QuadrantProblem()

        mesh, coords = Meshes.make_mesh((0.0, 1.0, 0.0, 1.0), nx, ny;
            boundaries=(:outflow, :outflow, :outflow, :outflow))
        Meshes.mark_contact_faces!(mesh, coords, (a, b) ->
            (abs(a[1]-0.5) < 1e-12 && abs(b[1]-0.5) < 1e-12) ||
            (abs(a[2]-0.5) < 1e-12 && abs(b[2]-0.5) < 1e-12))
        options = Options.DGOptions(degree=degree, OF=OF, characteristic=characteristic,
            positivity=positivity, cfl=cfl, max_retries = max_retries)

        return Discretization.DGCase(problem, mesh, coords, options; final_time=final_time)
    end

end # module QuadrantCase