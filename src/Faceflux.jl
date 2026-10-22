module FaceFlux

    using LinearAlgebra: I, dot

    import ..Euler
    import ..Geometry
    import ..Problems

    function outside_cons(inside_cons_vec, point_coords, boundary, normal, problem)
        if boundary == :wall
            return Euler.reflected_cons(inside_cons_vec, normal)
        elseif boundary == :fixed
            return Problems.initial_cons(problem, point_coords)
        elseif boundary == :outflow
            return copy(inside_cons_vec)
        end
        error("Unknown boundary: $(boundary)")
    end

end #