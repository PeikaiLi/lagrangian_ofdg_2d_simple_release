# 2D Lagrangian-Type DG/OFDG Solver 

An implementation of a high-order, cell-centered Lagrangian-type discontinuous Galerkin method for the 2D compressible Euler equations with **oscillation-free** option, based on the work of [Feng, Cheng, and Shu (2025)](https://link.springer.com/article/10.1007/s10915-025-03107-5).

Note: I just do an independent implementation so it is not the official code associated with the publication. We want to see whether this method works for some complex 2D cases.

The paper using  quadrilateral meshes so I use it as well to replicate and compare the result, but now I think a triangular meshes might be better.

## Features (so far)
Spatial Discretization
- Cell-centered DG formulation on moving quadrilateral meshes.
- Any polynomial approximation degree(I do degree in 0-3 which can be easily modify to any degree).
- HLLC-based numerical fluxes 
  <span style="color: red; background-color:white;">
  maybe one issue here dealing if right and left different gamma
  </span>
- Moving-mesh treatment of the DG basis functions.
- Special treatment of cell-average fluxes at marked interfaces and detected shocks.
- Pressure-jump-based shock detection.

Time Integration
+ Forward Euler.
+ Third-order strong stability preserving Runge–Kutta (SSPRK3).
+ fourth-order Runge–Kutta (RK4).
+ CFL-based adaptive time-step selection.
+ adaptive time-step selection if non positivity happen

Oscillation-Free
- Damping coefficients based on jumps of the DG solution and its derivatives.
- Optional local characteristic decomposition.
- Operator splitting between DG evolution and numerical damping.
- Exact integration of the discrete damping equation using a matrix exponential.
- Damping applied after each Runge–Kutta stage.

Postprocess
- Conservation test and convergence test for smooth problems and plotting converge rate.
- Visualization of moving meshes and numerical solutions (plot average value per cell so far, but in next version I will plot the detail in each cell better.)
- GIF for time evolution.

##  files
<!-- tree /F -->
+ `example/` is all the cases can be run in jupyter notebook. 
+ `design/ `is my notes about how I design my data structure and implement the method

