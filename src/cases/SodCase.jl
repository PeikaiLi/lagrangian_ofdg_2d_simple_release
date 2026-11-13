module SodCase

    import ..Euler
    import ..Problems
    import ..Meshes
    import ..Options
    import ..Discretization

    struct SodProblem <: Problems.EulerProblem
        name::Symbol
        smooth::Bool
    end

    SodProblem() = SodProblem(:sod, false)


    function Problems.initial_gamma(problem::SodProblem, point_coord) 
        x, y = point_coord

        # if (0-tol) <= x < 0.5
        if x < 0.5
            return 1.4
        end
        # if 0.5 <= x <= (1+tol)
        if 0.5 <= x
            return 5/3
        end

        throw(ArgumentError("esti gamma SodProblem requires 0 <= x <= 1; current x = $(x)"))
    end    


    function Problems.initial_cons(problem::SodProblem, point_coord)

        x, y = point_coord

        if  x < 0.5
            rho = 1.0
            u = 0.0
            v = 0.0
            p = 1.0
            return Euler.prim2cons(rho, u, v, p, 1.4)
        end

        if 0.5 <= x
            rho = 0.125
            u = 0.0
            v = 0.0
            p = 0.1
            return Euler.prim2cons(rho, u, v, p, 5/3)
        end

        throw(ArgumentError("esti init SodProblem requires 0 <= x <= 1; current x = $(x)"))
    end

    # shu use nx = 200 ny = 2;
    # then he claimed that he use 50*2
    # final_time = 0.2   
    # cfl unknow 
    function make_case(; nx=50, ny=2, degree=0, final_time=0.2,
                   OF=false, characteristic=false, positivity=false, cfl=0.3)

        iseven(nx) || throw(ArgumentError("Use even nx so x = 10 is a mesh face"))

        problem = SodProblem()
        # :wall is reflective boundary
        mesh, coords = Meshes.make_mesh((0.0, 1.0, 0.0, 0.1), nx, ny;
            boundaries=(:fixed, :fixed, :wall, :wall))

        Meshes.mark_contact_faces!(mesh, coords,
            (a, b) -> abs(a[1] - 0.5) < 1e-12 && abs(b[1] - 0.5) < 1e-12) #1E-12 might not be the best option?

        options = Options.DGOptions(degree=degree, OF=OF, characteristic=characteristic,
            positivity=positivity, cfl=cfl)

        return Discretization.DGCase(problem, mesh, coords, options; final_time=final_time)
    end

end # module SodCase