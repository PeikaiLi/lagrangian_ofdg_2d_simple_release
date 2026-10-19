module Problems

abstract type EulerProblem end

struct DGCase{P<:EulerProblem}
    initial::States.DGState
    mesh::Meshes.Mesh
    problem::P
    options::Options.DGOptions
end

end # module Problems