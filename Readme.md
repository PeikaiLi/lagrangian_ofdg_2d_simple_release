# 2D ALE DG Replication

So far, I've implemented the mesh vertex velocity calculation. Next, I'll refactor the integration routines to make the code clearer and more efficient.
After than I can finially update the moment and coefficient, then to run a complete sooth order test... I can't believe I was so Fucking slow.
Current files:

```text
.
├── README.md
└── code/
    ├── Options.jl          # Numerical settings
    ├── Meshes.jl           # Mesh connectivity, faces, and periodic pairs
    ├── Quadrature.jl       # Gauss rules and integration over quadrilaterals
    ├── Basis.jl            # Polynomial basis functions and their derivatives
    ├── Geometry.jl         # Cell geometry, mass matrices, and polynomial evaluation
    ├── States.jl           # DG state storage and initial projection
    ├── Euler.jl            # Euler variables, fluxes, and HLLC calculations
    ├── Problems.jl         # Interfaces for initial data, gamma, and source terms
    ├── Discretization.jl   # Builds a case from a problem, mesh, and options
    ├── FaceFlux.jl         # Face states, boundary conditions, and moving-face fluxes
    ├── MeshMotion.jl       # Vertex velocities and basis time derivatives
    └── cases/
        └── ConstantCase.jl # Uniform-flow example
```

Start with `make_case` in `code/cases/ConstantCase.jl` to see how the pieces fit together.

