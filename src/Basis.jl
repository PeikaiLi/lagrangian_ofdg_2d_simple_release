module Basis

function basis_powers(degree)
    exponents_ls = [(0, 0)]
    if degree >= 1
        append!(exponents_ls, [(1, 0), (0, 1)])
    end
    if degree >= 2
        append!(exponents_ls, [(2, 0), (0, 2), (1, 1)])
    end
    if degree >= 3
        append!(exponents_ls, [(3, 0), (0, 3), (2, 1), (1, 2)])
    end
    if degree < 0 || degree > 3
        error("degree must be 0, 1, 2 or 3")
    end
    return exponents_ls
end

function monomial_factor(x_power, y_power)
    # just match the paper, I actually don't think it matters at all.
    if (x_power == 2 && y_power == 0) || (x_power == 0 && y_power == 2)
        return 0.5
    end
    return 1.0
end

# We will do orthogonality phi0 = 1 later as paper required. For now, we will just do the monomials.
function raw_basis(point_coordinates, centroid_coordinates, hx, hy, exponents_ls; dx = 0, dy = 0)
    point = point_coordinates
    centroid = centroid_coordinates
    xi = 2.0 * (point[1] - centroid[1]) / hx
    eta = 2.0 * (point[2] - centroid[2]) / hy

    values = zeros(length(exponents_ls))

    for (k, (p, q)) in enumerate(exponents_ls)
        if p >= dx && q >= dy
            derivative_factor =
                factorial(p) / factorial(p - dx) *
                factorial(q) / factorial(q - dy)

            values[k] =
                monomial_factor(p, q) * derivative_factor *
                xi^(p - dx) * eta^(q - dy) *
                (2.0 / hx)^dx * (2.0 / hy)^dy
        end
    end

    return values
end

end # module Basis

