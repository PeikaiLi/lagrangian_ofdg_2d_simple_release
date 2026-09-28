module Basis

function basis_powers(degree)
    basis_power_pairs_list = [(0, 0)]
    if degree >= 1
        append!(basis_power_pairs_list, [(1, 0), (0, 1)])
    end
    if degree >= 2
        append!(basis_power_pairs_list, [(2, 0), (0, 2), (1, 1)])
    end
    if degree >= 3
        append!(basis_power_pairs_list, [(3, 0), (0, 3), (2, 1), (1, 2)])
    end
    if degree < 0 || degree > 3
        error("degree must be 0, 1, 2 or 3")
    end
    return basis_power_pairs_list
end

function monomial_factor(x_power, y_power)
    if (x_power == 2 && y_power == 0) || (x_power == 0 && y_power == 2)
        return 0.5
    end
    return 1.0
end

function raw_basis(point_coordinates, centroid_coordinates, hx, hy, basis_power_pairs_list; dx = 0, dy = 0)

    xi = 2.0 * (point_coordinates[1] - centroid_coordinates[1]) / hx
    eta = 2.0 * (point_coordinates[2] - centroid_coordinates[2]) / hy

    basis_values_vec = zeros(length(basis_power_pairs_list))

    for basis_index in eachindex(basis_power_pairs_list)
        x_power, y_power = basis_power_pairs_list[basis_index]
        if x_power >= dx && y_power >= dy
            derivative_factor = factorial(x_power) / factorial(x_power - dx) *
                factorial(y_power) / factorial(y_power - dy)

            basis_values_vec[basis_index] = monomial_factor(x_power, y_power) *
                derivative_factor * xi^(x_power - dx) * eta^(y_power - dy) * (2.0 / hx)^dx *
                (2.0 / hy)^dy
        end
    end

    return basis_values_vec
end

end # module Basis
