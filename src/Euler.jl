# this is similar to the 1D case.
module Euler
    import ..Options


    function prim2cons(rho, velocity_x, velocity_y, pressure_value, gamma)
        total_energy = pressure_value / (gamma - 1.0) + 0.5 * rho * (velocity_x^2 + velocity_y^2)

        cons_vec = [
            rho,
            rho * velocity_x,
            rho * velocity_y,
            total_energy,
        ]
        return cons_vec
    end

    function pressure(cons_vec, gamma)
        rho = cons_vec[1]
        kinetic_energy = (cons_vec[2]^2 + cons_vec[3]^2) / (2.0 * rho)
        return (gamma - 1.0) * (cons_vec[4] - kinetic_energy)
    end

    function cons2prim(cons_vec, gamma)
        rho = cons_vec[1]
        Options.check_positive(rho, "density must be positive")

        velocity_x = cons_vec[2] / rho
        velocity_y = cons_vec[3] / rho
        pressure_value = pressure(cons_vec, gamma)
        Options.check_positive(pressure_value, "pressure must be positive")

        return rho, velocity_x, velocity_y, pressure_value
    end

    function euler_flux(cons_vec, gamma)
        rho, velocity_x, velocity_y, pressure_value = cons2prim(cons_vec, gamma)

        # flux in x and y directions
        flux_x_vec = [
            rho * velocity_x,
            rho * velocity_x^2 + pressure_value,
            rho * velocity_x * velocity_y,
            (cons_vec[4] + pressure_value) * velocity_x,
        ]

        flux_y_vec = [
            rho * velocity_y,
            rho * velocity_x * velocity_y,
            rho * velocity_y^2 + pressure_value,
            (cons_vec[4] + pressure_value) * velocity_y,
        ]

        return flux_x_vec, flux_y_vec
    end

    # Riemann solver for the Euler equations



end # module Euler