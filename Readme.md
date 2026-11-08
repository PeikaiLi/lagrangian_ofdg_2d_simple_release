# 2D ALE DG Replication

So far,  Smooth test is done, which has attained the order as we expected. Next step will be "OF".
Current files:
<!-- tree /F -->
```text
C:\Users\Peikai_Li\Desktop\lagrangian_ofdg_2d_simple_release>tree /F
Folder PATH listing
Volume serial number is 2850-F372
C:.
│  .gitignore
│  Readme.md
│
├─design
│      n01meshes.pdf
│      n02_quadrature_mesh.pdf
│      n03Options&Basis&Cellgeometry.pdf
│      n04States .pdf
│      n05Euler&RiemannSolver.pdf
│      n06Problem&DGcase.pdf
│      n07caseIsentropic.pdf
│      n08meshmotion&basismotion.pdf
│      n09DGrhs.pdf
│
├─example
│      constant.ipynb
│      isentropic_case.ipynb
│
├─reference
│      19870013797_SSRK3.pdf
│      Gottlieb-Shu-Tadmor.SIREV-01_SSPRK.pdf
│      Shu_LagrangianTypeDG .pdf
│      sspreview.pdf
│      Yee_forRoe.pdf
│
├─results
│      constant_density.gif
│      constant_pressure.gif
│      constant_velocity_x.gif
│      constant_velocity_y.gif
│      IsentropicCasedensity.gif
│      IsentropicCasepressure.gif
│      IsentropicCasevelocity_x.gif
│      IsentropicCasevelocity_y.gif
│
├─src
│  │  Basis.jl
│  │  DGOperator.jl
│  │  Diagnostics.jl
│  │  Discretization.jl
│  │  Euler.jl
│  │  Faceflux.jl
│  │  Geometry.jl
│  │  load_all.jl
│  │  Meshes.jl
│  │  MeshMotion.jl
│  │  Options.jl
│  │  Positivity.jl
│  │  Problems.jl
│  │  Quadrature.jl
│  │  Solver.jl
│  │  States.jl
│  │  test_smooth_case.ipynb
│  │  TimeIntegration.jl
│  │  Visualization.jl
│  │
│  └─cases
│          ConstantCase.jl
│          IsentropicCase.jl
│
└─test_
        plot_test.ipynb
```

## It seems difficult and indeed it is difficult. But I find the best way is that do not only think how to code but just start it, starting to write on the paper and coding on computer. This can put my **RAM** clean and give me more space to think about next step.
