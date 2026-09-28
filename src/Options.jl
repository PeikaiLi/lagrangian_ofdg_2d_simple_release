module Options

struct DGOptions
    degree::Int
    use_ofdg::Bool
    use_characteristic_decomposition::Bool
    use_positivity_limiter::Bool
    cfl::Float64
    quadrature_order::Int
    density_floor::Float64
    pressure_floor::Float64
    max_retries::Int
    max_steps::Int
end

function DGOptions(; degree = 2, OF = true, characteristic = true,
                    positivity = false, cfl = nothing, quadrature_order = degree+2)
    if cfl === nothing
        if degree <= 1
            cfl = 0.5
        elseif degree == 2
            cfl = 0.3
        else
            cfl = 0.1
        end
    end
    return DGOptions(degree, OF, characteristic, positivity, Float64(cfl),
                     max(2, quadrature_order), 1.0e-13, 1.0e-13, 10, 50000)
end

function check_positive(value, message)
    if !isfinite(value) || value <= 0.0
        throw(DomainError(value, message))
    end
end

end # module Options
