# this is similar to the 1D case.
module Euler
    using LinearAlgebra: norm
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
    function hllc_star(cons_in,cons_out,gamma_in,gamma_out,out_normal)
        normal = out_normal / norm(out_normal)  # Normalize the normal vector
        tangent = [-normal[2], normal[1]]  # Tangent vector perpendicular to the normal
        # Compute the prim 
        rho_in, velocity_x_in, velocity_y_in, pressure_in = cons2prim(cons_in, gamma_in)
        rho_out, velocity_x_out, velocity_y_out, pressure_out = cons2prim(cons_out, gamma_out)

        # Compute the normal velocities
        velocity_n_in = velocity_x_in * normal[1] + velocity_y_in * normal[2]
        velocity_n_out = velocity_x_out * normal[1] + velocity_y_out * normal[2]

        velocity_t_in = velocity_x_in * tangent[1] + velocity_y_in * tangent[2]
        velocity_t_out = velocity_x_out * tangent[1] + velocity_y_out * tangent[2]

        # Compute the sound speeds
        c_in = sqrt(gamma_in * pressure_in / rho_in)
        c_out = sqrt(gamma_out * pressure_out / rho_out)

        # Roe averages
        D = sqrt(rho_out/rho_in)
        u_n_Roe = (velocity_n_in + D * velocity_n_out) / (1 + D)
        u_t_Roe = (velocity_t_in + D * velocity_t_out) / (1 + D)

        H_in = (cons_in[4] + pressure_in) / rho_in
        H_out = (cons_out[4] + pressure_out) / rho_out
        H_Roe = (H_in + D * H_out) / (1 + D)

        # I don't know if this is the best way. But I just put one here.
        gamma = (gamma_in + gamma_out) / 2 
        c_Roe = sqrt((gamma - 1) * (H_Roe - 0.5 * (u_n_Roe^2+u_t_Roe^2)))

        mu_mins = min(velocity_n_in - c_in, u_n_Roe - c_Roe)
        mu_plus = max(velocity_n_out + c_out, u_n_Roe + c_Roe)

        numerator_mu_star = 
            rho_out * velocity_n_out * (mu_plus - velocity_n_out) -
            rho_in * velocity_n_in * (mu_mins - velocity_n_in) +
            pressure_in - pressure_out
        denominator_mu_star = rho_out * (mu_plus - velocity_n_out) - rho_in * (mu_mins - velocity_n_in)
        mu_star = numerator_mu_star / denominator_mu_star

        p_star = rho_in * (velocity_n_in - mu_mins) * (velocity_n_in - mu_star) + pressure_in

        return mu_star, p_star
    end

end # module Euler