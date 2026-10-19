# this is similar to the 1D case.
module Euler
    import ..Options
    

    function is_admissible(cons_vec, gamma, options::Options.DGOptions)
        if !all(isfinite, cons_vec) ||
           cons_vec[1] <= options.density_floor
            return false

        pressure_value = pressure(cons_vec, gamma)
        if isfinite(pressure_value) 
            return false
        elseif pressure_value < options.pressure_floor
            return false
        end
        return true
    end

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

        normal = out_normal
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
        c_Roe_square = (gamma - 1) * (H_Roe - 0.5 * (u_n_Roe^2 + u_t_Roe^2))
        Options.check_positive(c_Roe_square, "HLLC Roe sound speed squared")
        c_Roe = sqrt(c_Roe_square)

        mu_mins = min(velocity_n_in - c_in, u_n_Roe - c_Roe)
        mu_plus = max(velocity_n_out + c_out, u_n_Roe + c_Roe)

        numerator_mu_star = 
            rho_out * velocity_n_out * (mu_plus - velocity_n_out) -
            rho_in * velocity_n_in * (mu_mins - velocity_n_in) +
            pressure_in - pressure_out
        denominator_mu_star = rho_out * (mu_plus - velocity_n_out) - rho_in * (mu_mins - velocity_n_in)
        Options.check_positive(denominator_mu_star, "HLLC zero denominator_mu_star in the moving basis")
        mu_star = numerator_mu_star / denominator_mu_star

        p_star = rho_in * (velocity_n_in - mu_mins) * (velocity_n_in - mu_star) + pressure_in

        return mu_star, p_star
    end

    # if boundary is a wall
    #, then the reflected cons_vec is used as a image state for the Riemann solver.
    # use HLLC to compute the flux.
    function reflected_cons(cons_vec, normal)
        reflect = copy(cons_vec)
        U = cons_vec[2:3]
        U_reflected = U - 2 * (U' * normal) * normal
        reflect[2:3] = U_reflected
            return reflect
    end

    # For OF we will do characteristic decomposition
    # I do it here
    function left_eigen_vector_Jacobians(ave_cons_vec, gamma, normal)
        rho, velocity_x, velocity_y, pressure = cons2prim(ave_cons_vec, gamma)
        c = sqrt(gamma * pressure / rho)

        u_tilde = velocity_x * normal[1] + velocity_y * normal[2]
        v_tilde = -velocity_x * normal[2] + velocity_y * normal[1]
       
        B_r = (gamma - 1) * ( 1/2 * (velocity_x^2 + velocity_y^2)) 

        R_inv = zeros(Float64, 4, 4)

        R_inv[1, 1] = 1/2 * (B_r + u_tilde*c)
        R_inv[2, 1] = v_tilde * c
        R_inv[3, 1] = c^2 - B_r
        R_inv[4, 1] = 1/2 * (B_r - u_tilde*c)

        R_inv[1, 2] = -1/2 *((gamma - 1)*velocity_x + normal[1]*c)
        R_inv[2, 2] = normal[2] * c
        R_inv[3, 2] = (gamma - 1) * velocity_x
        R_inv[4, 2] = -1/2 *((gamma - 1)*velocity_x - normal[1]*c) 

        R_inv[1, 3] = -1/2 *((gamma - 1)*velocity_y + normal[2]*c)
        R_inv[2, 3] = -normal[1] * c
        R_inv[3, 3] = (gamma - 1) * velocity_y
        R_inv[4, 3] = -1/2 *((gamma - 1)*velocity_y - normal[2]*c)

        R_inv[1, 4] = 1/2 *(gamma - 1)
        R_inv[2, 4] = 0
        R_inv[3, 4] = -(gamma - 1)
        R_inv[4, 4] = 1/2 *(gamma - 1)

        scaler = (gamma - 1) / c

        R_inv = scaler * R_inv
        return R_inv
    end

end # module Euler