# 2D ALE DG Replication

So far, I've implemented the mesh vertex velocity calculation. Next, I'll refactor the integration routines to make the code clearer and more efficient.
After than I can finially update the moment and coefficient, then to run a complete sooth order test... I can't believe I was so Fucking slow.
Current files:
<!-- tree /F -->
```text
C:.
│  Readme.md
│  
├─design
│      n01meshes.pdf
│      n02_quadrature_mesh.pdf
│      n03Options&Basis&Cellgeometry.pdf
│      n04States .pdf
│      n05Euler&RiemannSolver.pdf
│      n06Problem&DGcase.pdf
│      n08meshmotion&basismotion.pdf
│      
├─reference
│      Shu_LagrangianTypeDG .pdf
│      Yee_forRoe.pdf
│      
└─src
    │  Basis.jl
    │  DGOperator.jl
    │  Discretization.jl
    │  Euler.jl
    │  Faceflux.jl
    │  Geometry.jl
    │  Meshes.jl
    │  MeshMotion.jl
    │  Options.jl
    │  Problem.jl
    │  Quadrature.jl
    │  States.jl
    │  test.ipynb
    │  
    ├─.ipynb_checkpoints
    │      Untitled-checkpoint.ipynb
    │      
    └─cases
            ConstantCase.jl
```

Start with `make_case` in `code/cases/ConstantCase.jl` 



## It seems difficult and indeed it is difficult. But I find the best way is that do not only think how to code but just start it, starting to write on the paper and coding on computer. This can put my **RAM** clean and give me more space to think about next step.
