module LaxCase

    import ..Euler
    import ..Problems
    import ..Meshes
    import ..Options
    import ..Discretization

    struct LaxProblem <: Problems.EulerProblem
        name::Symbol
        smooth::Bool
    end

    LaxProblem() = LaxProblem(:lax, false)

    Problems.initial_gamma(problem::LaxProblem, point_coord) = 1.4

    function Problems.initial_cons(problem::LaxProblem, point_coord)
        x, y = point_coord
        if 5 <= x < 10
            rho = 0.445
            u = 0.698
            v = 0.0
            p = 3.528
        else
            rho = 0.5
            u = 0.0
            v = 0.0
            p = 0.571
        end
        return Euler.prim2cons(rho, u, v, p, 1.4)
    end

    # shu use nx = 100
    # final_time = 1    
    # cfl unknow 
    function make_case(; nx=20, ny=2, degree=2, final_time=0.2,
                   OF=false, characteristic=false, positivity=false, cfl=0.3)

        iseven(nx) || throw(ArgumentError("Use even nx so x = 10 is a mesh face"))

        problem = LaxProblem()
        # :wall is reflective boundary
        mesh, coords = Meshes.make_mesh((5, 15, 0.0, 0.1), nx, ny;
            boundaries=(:fixed, :fixed, :wall, :wall))

        Meshes.mark_contact_faces!(mesh, coords,
            (a, b) -> abs(a[1] - 10) < 1e-12 && abs(b[1] - 10) < 1e-12)

        options = Options.DGOptions(degree=degree, OF=OF, characteristic=characteristic,
            positivity=positivity, cfl=cfl)

        return Discretization.DGCase(problem, mesh, coords, options; final_time=final_time)
    end

end # module LaxCase