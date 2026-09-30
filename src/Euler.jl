module Euler

import ..Options

# Euler equations in 2D
function prim_to_cons(rho, velocity_x, velocity_y, pressure_value, gamma)
    total_energy = pressure_value / (gamma - 1.0) + 0.5 * rho * (velocity_x^2 + velocity_y^2)

    state_vec = [
        rho,
        rho * velocity_x,
        rho * velocity_y,
        total_energy,
    ]
    return state_vec
end

function pressure(state_vec, gamma)
    rho = state_vec[1]
    kinetic_energy = (state_vec[2]^2 + state_vec[3]^2) / (2.0 * rho)
    return (gamma - 1.0) * (state_vec[4] - kinetic_energy)
end

function cons_to_prim(state_vec, gamma)
    rho = state_vec[1]
    Options.check_positive(rho, "density must be positive")

    velocity_x = state_vec[2] / rho
    velocity_y = state_vec[3] / rho
    pressure_value = pressure(state_vec, gamma)
    Options.check_positive(pressure_value, "pressure must be positive")

    return rho, velocity_x, velocity_y, pressure_value
end

function euler_flux(state_vec, gamma)
    rho, velocity_x, velocity_y, pressure_value = cons_to_prim(state_vec, gamma)

    flux_x_list = [
        rho * velocity_x,
        rho * velocity_x^2 + pressure_value,
        rho * velocity_x * velocity_y,
        (state_vec[4] + pressure_value) * velocity_x,
    ]

    flux_y_list = [
        rho * velocity_y,
        rho * velocity_x * velocity_y,
        rho * velocity_y^2 + pressure_value,
        (state_vec[4] + pressure_value) * velocity_y,
    ]

    return flux_x_list, flux_y_list
end

function hllc_star(left_state, right_state, left_gamma, right_gamma, normal)

    rho_left, velocity_x_left, velocity_y_left, pressure_left = cons_to_prim(left_state, left_gamma)
    rho_right, velocity_x_right, velocity_y_right, pressure_right =
        cons_to_prim(right_state, right_gamma)

    normal_velocity_left = velocity_x_left * normal[1] + velocity_y_left * normal[2]
    normal_velocity_right = velocity_x_right * normal[1] + velocity_y_right * normal[2]

    sound_speed_left = sqrt(left_gamma * pressure_left / rho_left)
    sound_speed_right = sqrt(right_gamma * pressure_right / rho_right)

    roe_weight_left = sqrt(rho_left)
    roe_weight_right = sqrt(rho_right)
    roe_weight_sum = roe_weight_left + roe_weight_right

    roe_normal_velocity = (roe_weight_left * normal_velocity_left +
         roe_weight_right * normal_velocity_right) / roe_weight_sum
    roe_sound_speed = (roe_weight_left * sound_speed_left +
         roe_weight_right * sound_speed_right) / roe_weight_sum

    wave_speed_minus = min(
        normal_velocity_left - sound_speed_left,
        roe_normal_velocity - roe_sound_speed,
    )
    wave_speed_plus = max(
        normal_velocity_right + sound_speed_right,
        roe_normal_velocity + roe_sound_speed,
    )

    denominator = rho_right * (wave_speed_plus - normal_velocity_right) -
        rho_left * (wave_speed_minus - normal_velocity_left)
    Options.check_positive(denominator, "HLLC denominator")

    star_normal_velocity = (rho_right * normal_velocity_right *
         (wave_speed_plus - normal_velocity_right) - rho_left * normal_velocity_left *
         (wave_speed_minus - normal_velocity_left) + pressure_left - pressure_right) / denominator

    star_pressure = pressure_left + rho_left * (normal_velocity_left - wave_speed_minus) *
        (normal_velocity_left - star_normal_velocity)

    return star_normal_velocity, star_pressure
end



function reflected_state(inside_state, normal, wall_normal_speed)
    normal_momentum = inside_state[2] * normal[1] + inside_state[3] * normal[2]

    correction = 2.0 * (inside_state[1] * wall_normal_speed - normal_momentum)

    rho = inside_state[1]
    momentum_x = inside_state[2] + correction * normal[1]
    momentum_y = inside_state[3] + correction * normal[2]
    total_energy = inside_state[4] +
        2.0 * inside_state[1] * wall_normal_speed^2 - 2.0 * wall_normal_speed * normal_momentum

    return [rho, momentum_x, momentum_y, total_energy]
end

function wave_curve(trial_pressure, rho, reference_pressure, gamma)
    # Compute the wave curve for a given trial pressure, reference density, reference pressure, and specific heat ratio (gamma).
    sound_speed = sqrt(gamma * reference_pressure / rho)
    if trial_pressure > reference_pressure
        coefficient_a = 2.0 / ((gamma + 1.0) * rho)
        coefficient_b = (gamma - 1.0) / (gamma + 1.0) * reference_pressure
        return (trial_pressure - reference_pressure) *
               sqrt(coefficient_a / (trial_pressure + coefficient_b))
    end

    exponent = (gamma - 1.0) / (2.0 * gamma)
    return 2.0 * sound_speed / (gamma - 1.0) * ((trial_pressure / reference_pressure)^exponent - 1.0)
end

function left_characteristic_matrix(state_vec, gamma, normal)
    rho, velocity_x, velocity_y, pressure_value = cons_to_prim(state_vec, gamma)

    sound_speed = sqrt(gamma * pressure_value / rho)
    gamma_minus_one = gamma - 1.0
    normal_velocity = velocity_x * normal[1] + velocity_y * normal[2]
    tangent_velocity = -velocity_x * normal[2] + velocity_y * normal[1]
    kinetic_factor = 0.5 * gamma_minus_one * (velocity_x^2 + velocity_y^2)

    left_characteristic_mat = zeros(4, 4)
    left_characteristic_mat[1, :] = [
        0.5 * (kinetic_factor + normal_velocity * sound_speed),
        -0.5 * (gamma_minus_one * velocity_x + normal[1] * sound_speed),
        -0.5 * (gamma_minus_one * velocity_y + normal[2] * sound_speed),
        0.5 * gamma_minus_one,
    ]
    left_characteristic_mat[2, :] = [
        tangent_velocity * sound_speed,
        normal[2] * sound_speed,
        -normal[1] * sound_speed,
        0.0,
    ]
    left_characteristic_mat[3, :] = [
        sound_speed^2 - kinetic_factor,
        gamma_minus_one * velocity_x,
        gamma_minus_one * velocity_y,
        -gamma_minus_one,
    ]
    left_characteristic_mat[4, :] = [
        0.5 * (kinetic_factor - normal_velocity * sound_speed),
        -0.5 * (gamma_minus_one * velocity_x - normal[1] * sound_speed),
        -0.5 * (gamma_minus_one * velocity_y - normal[2] * sound_speed),
        0.5 * gamma_minus_one,
    ]

    return (gamma_minus_one / sound_speed) * left_characteristic_mat
end

function is_admissible(state_vec, gamma, options)
    if !all(isfinite, state_vec) ||
       state_vec[1] <= options.density_floor
        return false
    end
    return pressure(state_vec, gamma) > options.pressure_floor
end

end # module Euler
