module Options

    struct DGOptions
        degree::Int
        quadrature_order::Int # Degree of precision of the 1D Gauss rule.
        cfl::Float64

        use_ofdg::Bool
        use_characteristic_decomposition::Bool

        use_positivity_limiter::Bool
        density_floor::Float64
        pressure_floor::Float64
        

        max_retries::Int
        max_steps::Int

        tol::Float64
    end

    function DGOptions(; degree = 2, OF = true, characteristic = true,
                        positivity = false, cfl = nothing, quadrature_order = 2 * degree + 3)
        if cfl === nothing
            if degree <= 1
                cfl = 0.5
            elseif degree == 2
                cfl = 0.3
            else
                cfl = 0.1
            end
        end
        return DGOptions(
            degree,
            max(3, quadrature_order),
            Float64(cfl),
            OF,
            characteristic,
            positivity,
            1.0e-13,
            1.0e-13,
            10,
            50000,
            1.0e-13,
        )
    end

    function check_positive(value, message)
        if !isfinite(value) || value <= 0.0
            throw(DomainError(value, message))
        end
    end

end # module Options
