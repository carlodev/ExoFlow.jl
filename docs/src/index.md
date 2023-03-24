# ExoFlow.jl
*Documentation of ExoFlow.jl for solving incompressible Navier-Stokes using stabilized Finite Element Method, in specific Streamline-Upwind Petrov-Galerkin (SUPG) and Variational Multiscale Method (VMS)*

## Installation
At the moment it is not a registered Julia package. For installing, from the `REPL` just press `]`.
```example
(@1.8) pkg> add https://github.com/carlodev/ExoFlow.jl 
```

For a complete and smooth experience is suggested to install also the free software [ParaView](https://www.paraview.org/) which allows to graphically visualize the results.

## Introduction



## Package General Features
- Implementation of SUPG and VMS formulation for same order elements for velocity and pressure
- Use of different ode method: ``θ-method`` with different values of ``θ`` for velocity and pressure, and ``α-mehtod``
- Using custom Meshes created with [`gmsh`](https://gmsh.info/). For airfoils the package `AirfoilGmsh.jl` has been developed for speeding up the process
- Solve 2D and 3D cases
- Possibility of chosing the backend thanks to `PartitionedArrays.jl`. It can be run in the `REPL` or in `MPI`

## Known Issues
- The gamg preoconditioner is not compatible with `SequentialBackend()`
- The final results is a .vtu file, which can be open with ParaView. It is an index file of .pvtu files in `Results/` folder. Each .pvtu file is the solution at a single time step and each .pvtu is and index file of .pvtu files. Each of this .vtu files is the portion of one processor of the solution of a spefic time step. The .pvtu files have a wrong indexation, the user has to remove the `Results/`.

## Acknowledgement
The Variational Multiscale formualation is based on the paper: Bazilevs Y. et al., *Variational multiscale residual-based turbulence modeling forlarge eddy simulation of incompressible flows* , 10.1016/j.cma.2007.07.016

The Streamline Upwind Petrov Galerkin is based on the PhD thesis:  Tamas Banyai,  *Development of Stabilized Finite ElementMethod for Numerical Simulation ofTurbulent Incompressible Single andEulerian-Eulerian Two-Phase Flows* 

The linearization of the equations is based on the PhD thesis:  Bart Janssen, *Numerical modeling and experimentalinvestigation of fine particle coagulationand dispersion in dilute flows* 

The package is based on the use of [`Gridap`](https://github.com/gridap/Gridap.jl): Santiago Badia and Francesc Verdugo,  *Gridap: An extensible Finite Element toolbox in Julia*, 10.21105/joss.02520

The solver used in MPI is PETSc: Satish Balay et al., https://petsc.org/