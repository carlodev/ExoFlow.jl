# Package usage
The package allows the user to set a wide variety of options.
Problem Settings:
- `:N` = number of divisions for each dimension. The example creates a 100x100 grid.
- `:D` = dimension. It can be 2 or 3.
- `:order` = order of the elements. At the moment just the order 1 is tested.
- `:case` it can be `"TaylorGreen", "LidDriven", "Cylinder", "Channel", "Airfoil"`.
- `:u_in` the inlet velocity for `"Airfoil"` and `"Cylinder"`, or the lid velocity for `"LidDriven`
- `:c` chord length in the `"Airfoil"` case, or dimension of lid for `"LidDriven"`. It is used to compute the viscosity `:ν` from the Reynolds and velocity
- `:Re` Reynolds number. 
- `:ν` kinematic viscosity. It can be overwritten in order to satisfy the Reynolds, in this case a warning informs the user.
- `:ρ` density. It used just to compute the force. The advice is to keep it `1.0` and just set the Reynolds.
- `:body_force` is non-zero generally just for the case of a periodic channel.
- `:periodic` used only in the `"Channel"` case. It can be set `true` or `false`


Time settings
- `:t0` = starting time.
- `:dt` = time step length.
- `:tF` = end time.
- `:t_endramp` = for high reynolds cases, like airfoils and lid driven, for improving numeric stability the inlet velocity (or the lid velocity) are increased from 0 up the desired value in the time between :t0 and :t_endramp. If :t0 = :t_endramp there is no ramping.

Ode Settings
- `:ode_method` = can be `:ThetaMethod` or `:AlphaMethod`.
- ``:θ`` = parameter required if `:ThetaMethod` is selected. It can a be single float, es `1.0` or a vector `[0.5, 1.0]`. In the case in the example it means that a ``θ = 0.5`` is used for velocity and ``θ = 1.0`` is used for pressure.
- `:ρ∞` = parameter required if `:ThetaMethod` is selected, it controls the numeric viscosity. A suggested value is `ρ∞ = 0.8`

Numeric Settings
- `:solver` for solving the equations, it can be `:julia` or `:petsc`. The general advice is to use `:petsc` expecially in `MPI`
- `:method` can be `:SUPG` or `:VMS`
- `:Cᵢ` is a vector containing stabilization coefficients used for the `:VMS`. The suggested values are `[4,36]`, 10.1016/j.compfluid.2008.10.003
- `:options` the settings for the petsc solver. It call the function `petsc_options(args)`, where `args` can be `:ksplu`, `:kspgamg` for linear case, `:sneslu`, `:snesgamg` for a non linear case. For a more detail explanation [`petsc_options`](@ref)
- `:linear` can be `true` or `false`. It linearizes the convective term using a Taylor expansion
- `:steady` is set to `false`. At the moment is not implemented `steady` solution of the equations 

Print Settings
- `:printmodel` can be `true` or `false`. If `true` mesh is saved as a .pvtu file.
- `:printinitial` can be `true` or `false`. If `true` saves the flowfield at `t0`. It is useful when restarting from a previous solution.
- `:benchmark_mode` => can be `true` or `false`. If `true` it does not print the solution, and it gives the time needed for computing the iteration form the 2nd till the end. The first iteration is not taken into account for computing the time because of precompilation.
- `:print_last` can be `true`or `false`. If true the solution is printed only for the last time step, if false, solution files are created for each time step. It is useful when you want to check the final state (`LidDriven` or `TaylorGreen`).
 
Mesh Settings
- `:mesh_gen` can be `true` or `false`. It has to be set `true` for cases where the mesh is read from a `.msh` file.
- `:mesh_file` is a string with the name of the `.msh` that can be read. By default it points to the folder `/models` of the package.


Partitioning Settings
- `:np_i` set the number of division in the `i` axes for carthesian problems. For non cartesian problems it does not matter how the cores are split into the dimensions as long as:
    - in the 2D case `:np_x * :np_y` has to equal to the `MPI` ranks.
    - in the 3D case `:np_x * :np_y * :np_z` has to equal to the `MPI` ranks.
- `backend` can be `MPIBackend()` or `SequentialBackend()`. 

Restarting Settings
- `:restart` can be `true` or `false`. If `false` the initial conditions are computed internally using `:u_in` or analytical solution (`"TaylorGreen"`). 
- `:restart_file` is used only if `:restart`is true. It is a `.csv` file created from ParaView using the SpreadSheet. It has the list of fo velocity and pressure in each node. It is better to run `clean grid` in Paraview before for get rid of duplicate points.

Turbulence Settings
For creating turbulence the package [`SyntheticEddyMethod`](https://github.com/carlodev/SyntheticEddyMethod.jl) is used.
- `:start_condition` for the channel, still work in progress.
- `:TI` Turbulence Intensity for the inlet. If it is set `0.0` it means no turbulence.
- `:Vbox => turbulence_box()` contains the information of the virtual box where the Eddies are created. More details in the documentation of [`SyntheticEddyMethod`](https://github.com/carlodev/SyntheticEddyMethod.jl). The parameters can be adjsuted in `Turbulence_Settings.jl` file. 
- `:Re_filename` it contains the string of the Reynolds stress file, which is a `.xlsx` file. If you want to create turbulence from the `:TI` parameters set it to `"none"`


```Example
using ExoFlow
using PartitionedArrays
include("Turbulence_Settings.jl")

backend = SequentialBackend() #or MPIBackend() SequentialBackend()
function instantiate_parameters()
    parameters = Dict(
    :N => 100,
    :D => 2, #Dimension
    :order => 1,
    
    :t0 => 0.0,
    :dt => 0.01,
    :tF => 0.05,
    :t_endramp => 0.0,

    :case => "TaylorGreen",
    :solver => :petsc,
    :method => :VMS,
    :ode_method => :ThetaMethod,
    :θ=>[0.5, 1.0],
    :ρ∞ => 1.0,
    :Re => 1_000,
    :c => 1, #chord lenght [m], used for naca (=1), cylinder (=0.1), liddriven = (1), 0.5
    :u_in => 1.0,  # =1.0 for lid driven 
    
    :periodic => false,

    :printmodel => false,
    
    :mesh_gen => false,

    :linear => true,
    :steady => false,

    :debug_mode => false,
    :benchmark_mode => false,
    :printinitial => false,
    :print_last =>false,
    :mesh_file => "DU89_i_2D.msh",
    :Cᵢ => [4, 36],    
    :options => petsc_options(:ksplu),

    :nls_trace =>true,
    :nls_iter => 20,

    :ν => 1.0e-5,  #channel = 0.0001472, 
    :ρ => 1.0, #kg/m3 density
    :body_force => 0.0, #channel = 0.00337204
    
    :np_x => 2, #number of processors in X
    :np_y => 2, #number of processors in Y
    :np_z => 1, #number of processors in Z

    :start_condition => :turbulent,
    :restart => false,
    :restart_file => "DU89_0p85.csv",
    :TI =>0.01,
    :Vbox => turbulence_box(),
    :Re_filename => "none"
      )

    parameters = initialize_parameters(parameters)
    
    return parameters
end

params = instantiate_parameters()


ExoFlow.main((params, backend))
```

The example above run in the `REPL` emulating a parallel run over 4 processors (you can see it by the options `:np_x` and `:np_y`). 
The final results is a .pvd file, which can be open with ParaView. It is an index file of .pvtu files in `Results/` folder. Each .pvtu file is the solution at a single time step and each .pvtu is and index file of .pvtu files. Each of this .vtu files is the portion of one processor of the solution of a spefic time step.

!!! info "numeric info" 
    Changing the backend to `MPIBackend()` allows it to run in MPI. 
