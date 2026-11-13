module ShuOsherCase
# I choose the
# boundaries = (:fixed, :outflow, :wall, :wall), 
# I have no idea what shu use.
    import ..Euler
    import ..Problems
    import ..Meshes
    import ..Options
    import ..Discretization


    struct ShuOsherProblem <: Problems.EulerProblem
        name::Symbol
        smooth::Bool
    end

    ShuOsherProblem() = ShuOsherProblem(:shu_osher, false)


    function Problems.initial_gamma(
        problem::ShuOsherProblem,
        point_coord,
    )
        return 1.4
    end


    function Problems.initial_cons(
        problem::ShuOsherProblem,
        point_coord,
    )
        x = point_coord[1]

        if x < -4.0
            density = 3.857143
            velocity_x = 2.629369
            pressure_value = 10.33333
        else
            density = 1.0 + 0.2 * sin(5.0 * x)
            velocity_x = 0.0
            pressure_value = 1.0
        end

        return Euler.prim2cons(
            density,
            velocity_x,
            0.0,
            pressure_value,
            1.4,
        )
    end


    function make_case(;
        nx = 100,
        ny = 2,
        degree = 0,
        final_time = 1.8,
        OF = true,
        characteristic = true,
        positivity = false,
        cfl = 0.1,
        max_retries = 4,
    )
        if degree == 0
            max_retries = 0
        end
        
        # Align a cell face with the initial discontinuity at x = -4.
        if nx <= 0 || nx % 5 != 0
            throw(ArgumentError(
                "nx must be positive and divisible by 5",
            ))
        end

        problem = ShuOsherProblem()

        mesh, coords = Meshes.make_mesh(
            (-10.0, 5.0, 0.0, 0.1),
            nx,
            ny;
            boundaries = (:fixed, :outflow, :wall, :wall),
        )

        # Retain the initial-interface label, following the Sod/Lax cases.
        # This initial jump contains a shock, not a pure contact.
        Meshes.mark_contact_faces!(
            mesh,
            coords,
            (a, b) ->
                abs(a[1] + 4.0) < 1e-12 &&
                abs(b[1] + 4.0) < 1e-12,
        )

        options = Options.DGOptions(
            degree = degree,
            OF = OF,
            characteristic = characteristic,
            positivity = positivity,
            cfl = cfl,
            max_retries = max_retries,
        )

        return Discretization.DGCase(
            problem,
            mesh,
            coords,
            options;
            final_time = final_time,
        )
    end

end # module ShuOsherCase