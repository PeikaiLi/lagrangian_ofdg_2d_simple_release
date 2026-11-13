module IsentropicCase

    import ..Euler
    import ..Problems
    import ..Meshes
    import ..Options
    import ..Discretization

    struct IsentropicProblem <: Problems.EulerProblem
        name::Symbol
        smooth::Bool
    end
    IsentropicProblem() = IsentropicProblem(:isentropic, true)

    Problems.initial_gamma(problem::IsentropicProblem, point_coord) = 3.0

    function Problems.initial_cons(problem::IsentropicProblem, point_coord)
        beta = 1.0 + 0.2 * sin((point_coord[1] + point_coord[2]) / 2.0)
        rho = beta / sqrt(6.0)
        return Euler.prim2cons(rho, beta / 2.0, beta / 2.0, rho^3, 3.0)
    end

function Problems.exact_cons(problem::IsentropicProblem, point_coord, time)
    (isfinite(time) && 0.0 <= time < 5.0) || throw(ArgumentError("The smooth exact solution requires 0 <= time < 5"))
    phase = (point_coord[1] + point_coord[2]) / 2.0
    # beta is in [0.8, 1.2]; the residual is strictly increasing before t=5.
    lower_beta, upper_beta = 0.8, 1.2
    beta = 1.0 + 0.2 * sin(phase)
    for iteration_index in 1:80
        angle = phase - beta * time
        residual = beta - 1.0 - 0.2 * sin(angle)
        if abs(residual) <= 1.0e-14
            rho = beta / sqrt(6.0)
            return Euler.prim2cons(rho, beta / 2.0, beta / 2.0, rho^3, 3.0)
        end
        if residual > 0.0
            upper_beta = beta
        else
            lower_beta = beta
        end
        derivative = 1.0 + 0.2 * time * cos(angle)
        trial_beta = beta - residual / derivative
        beta = lower_beta < trial_beta < upper_beta ? trial_beta : (lower_beta + upper_beta) / 2.0
    end
    error("The implicit exact solution did not converge")
end

function make_case(; nx=4, ny=nx, degree=2, OF=false, positivity=false, final_time=0.1, cfl=0.1)
    
    (isfinite(final_time) && 0.0 <= final_time < 5.0) || throw(ArgumentError("Need final_time < 5 for the smooth test"))

    problem = IsentropicProblem()

    mesh, vertices_coords_mat = Meshes.make_mesh(
        (0.0, 4.0*pi, 0.0, 4.0*pi), nx, ny;
        boundaries=(:periodic, :periodic, :periodic, :periodic),
    )
    
    options = Options.DGOptions(degree=degree, OF=OF, positivity=positivity, cfl=cfl)

    return Discretization.DGCase(problem, mesh, vertices_coords_mat, options; final_time=final_time)
end

end # module IsentropicCase
