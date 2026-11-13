module BlastWaveCase

    import ..Euler
    import ..Problems
    import ..Meshes
    import ..Options
    import ..Discretization


    struct BlastWaveProblem <: Problems.EulerProblem
        name::Symbol
        smooth::Bool
    end

    BlastWaveProblem() = BlastWaveProblem(:blast_wave, false)


    function Problems.initial_gamma(
        problem::BlastWaveProblem,
        point_coord,
    )
        return 1.4
    end


    function Problems.initial_cons(
        problem::BlastWaveProblem,
        point_coord,
    )
        x = point_coord[1]

        if x < 0.1
            pressure_value = 1000.0
        elseif x < 0.9
            pressure_value = 0.01
        else
            pressure_value = 100.0
        end

        return Euler.prim2cons(
            1.0,              # density
            0.0,              # velocity_x
            0.0,              # velocity_y
            pressure_value,
            1.4,
        )
    end


    function make_case(;
        nx = 100,
        ny = 2,
        degree = 0,
        final_time = 0.038,
        OF = true,
        characteristic = true,
        positivity = false,
        cfl = 0.1,
        max_retries = 4,
    )
        # Align cell faces with the two initial discontinuities.
        if nx <= 0 || nx % 10 != 0
            throw(ArgumentError(
                "nx must be positive and divisible by 10",
            ))
        end

        problem = BlastWaveProblem()

        mesh, coords = Meshes.make_mesh(
            (0.0, 1.0, 0.0, 0.01),
            nx,
            ny;
            boundaries = (:wall, :wall, :wall, :wall),
        )

        # Retain initial-interface labels, following the existing Sod/Lax cases.
        # The field is named contact, although these states have pressure jumps.
        Meshes.mark_contact_faces!(
            mesh,
            coords,
            (a, b) ->
                (
                    abs(a[1] - 0.1) < 1e-12 &&
                    abs(b[1] - 0.1) < 1e-12
                ) ||
                (
                    abs(a[1] - 0.9) < 1e-12 &&
                    abs(b[1] - 0.9) < 1e-12
                ),
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

end # module BlastWaveCase