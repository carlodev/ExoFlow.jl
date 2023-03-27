var documenterSearchIndex = {"docs":
[{"location":"mpi/#MPI-Run","page":"MPI Run","title":"MPI Run","text":"","category":"section"},{"location":"mpi/","page":"MPI Run","title":"MPI Run","text":"The code can be run in Message Passing Interface (MPI).  The code is made in such a way that it can run:","category":"page"},{"location":"mpi/","page":"MPI Run","title":"MPI Run","text":"in the REPL, selecting backend = SequentialBackend()\nin MPI, selecting backend = MPIBackend()","category":"page"},{"location":"mpi/","page":"MPI Run","title":"MPI Run","text":"When running in MPI the code cannot be easily executed in the REPL. Instead, one has to run them from a terminal using the mpiexecjl script as provided by MPI.jl. ","category":"page"},{"location":"mpi/","page":"MPI Run","title":"MPI Run","text":"From the folder /run for example with the command: ","category":"page"},{"location":"mpi/","page":"MPI Run","title":"MPI Run","text":"mpiexecjl -n 4 julia --project=../ run_ExoFlow.jl","category":"page"},{"location":"usage/#Package-usage","page":"Usage","title":"Package usage","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"The package allows the user to set a wide variety of options. Problem Settings:","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":N = number of divisions for each dimension. The example creates a 100x100 grid.\n:D = dimension. It can be 2 or 3.\n:order = order of the elements. At the moment just the order 1 is tested.\n:case it can be \"TaylorGreen\", \"LidDriven\", \"Cylinder\", \"Channel\", \"Airfoil\".\n:u_in the inlet velocity for \"Airfoil\" and \"Cylinder\", or the lid velocity for \"LidDriven\n:c chord length in the \"Airfoil\" case, or dimension of lid for \"LidDriven\". It is used to compute the viscosity :ν from the Reynolds and velocity\n:Re Reynolds number. \n:ν kinematic viscosity. It can be overwritten in order to satisfy the Reynolds, in this case a warning informs the user.\n:ρ density. It used just to compute the force. The advice is to keep it 1.0 and just set the Reynolds.\n:body_force is non-zero generally just for the case of a periodic channel.\n:periodic used only in the \"Channel\" case. It can be set true or false","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Time settings","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":t0 = starting time.\n:dt = time step length.\n:tF = end time.\n:t_endramp = for high reynolds cases, like airfoils and lid driven, for improving numeric stability the inlet velocity (or the lid velocity) are increased from 0 up the desired value in the time between :t0 and :tendramp. If :t0 = :tendramp there is no ramping.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Ode Settings","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":ode_method = can be :ThetaMethod or :AlphaMethod.\nθ = parameter required if :ThetaMethod is selected. It can a be single float, es 1.0 or a vector [0.5, 1.0]. In the case in the example it means that a θ = 05 is used for velocity and θ = 10 is used for pressure.\n:ρ∞ = parameter required if :ThetaMethod is selected, it controls the numeric viscosity. A suggested value is ρ∞ = 0.8","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Numeric Settings","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":solver for solving the equations, it can be :julia or :petsc. The general advice is to use :petsc expecially in MPI\n:method can be :SUPG or :VMS\n:Cᵢ is a vector containing stabilization coefficients used for the :VMS. The suggested values are [4,36], 10.1016/j.compfluid.2008.10.003\n:options the settings for the petsc solver. It call the function petsc_options(args), where args can be :ksplu, :kspgamg for linear case, :sneslu, :snesgamg for a non linear case. For a more detail explanation petsc_options\n:linear can be true or false. It linearizes the convective term using a Taylor expansion\n:steady is set to false. At the moment is not implemented steady solution of the equations ","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Print Settings","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":printmodel can be true or false. If true mesh is saved as a .pvtu file.\n:printinitial can be true or false. If true saves the flowfield at t0. It is useful when restarting from a previous solution.\n:benchmark_mode => can be true or false. If true it does not print the solution, and it gives the time needed for computing the iteration form the 2nd till the end. The first iteration is not taken into account for computing the time because of precompilation.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Mesh Settings","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":mesh_gen can be true or false. It has to be set true for cases where the mesh is read from a .msh file.\n:mesh_file is a string with the name of the .msh that can be read. By default it points to the folder /models of the package.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Partitioning Settings","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":np_i set the number of division in the i axes for carthesian problems. For non cartesian problems it does not matter how the cores are split into the dimensions as long as:\nin the 2D case :np_x * :np_y has to equal to the MPI ranks.\nin the 3D case :np_x * :np_y * :np_z has to equal to the MPI ranks.\nbackend can be MPIBackend() or SequentialBackend(). ","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Restarting Settings","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":restart can be true or false. If false the initial conditions are computed internally using :u_in or analytical solution (\"TaylorGreen\"). \n:restart_file is used only if :restartis true. It is a .csv file created from ParaView using the SpreadSheet. It has the list of fo velocity and pressure in each node. It is better to run clean grid in Paraview before for get rid of duplicate points.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Turbulence Settings For creating turbulence the package SyntheticEddyMethod is used.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":":start_condition for the channel, still work in progress.\n:TI Turbulence Intensity for the inlet. If it is set 0.0 it means no turbulence.\n:Vbox => turbulence_box() contains the information of the virtual box where the Eddies are created. More details in the documentation of SyntheticEddyMethod. The parameters can be adjsuted in Turbulence_Settings.jl file. \n:Re_filename it contains the string of the Reynolds stress file, which is a .xlsx file. If you want to create turbulence from the :TI parameters set it to \"none\"","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"using ExoFlow\nusing PartitionedArrays\ninclude(\"Turbulence_Settings.jl\")\n\nbackend = SequentialBackend() #or MPIBackend() SequentialBackend()\nfunction instantiate_parameters()\n    parameters = Dict(\n    :N => 100,\n    :D => 2, #Dimension\n    :order => 1,\n    \n    :t0 => 0.0,\n    :dt => 0.01,\n    :tF => 0.05,\n    :t_endramp => 0.0,\n\n    :case => \"TaylorGreen\",\n    :solver => :petsc,\n    :method => :VMS,\n    :ode_method => :ThetaMethod,\n    :θ=>[0.5, 1.0],\n    :ρ∞ => 1.0,\n    :Re => 1_000,\n    :c => 1, #chord lenght [m], used for naca (=1), cylinder (=0.1), liddriven = (1), 0.5\n    :u_in => 1.0,  # =1.0 for lid driven \n    \n    :periodic => false,\n\n    :printmodel => false,\n    \n    :mesh_gen => false,\n\n    :linear => true,\n    :steady => false,\n\n    :debug_mode => false,\n    :benchmark_mode => false,\n    :printinitial => false,\n\n    :mesh_file => \"DU89_i_2D.msh\",\n    :Cᵢ => [4, 36],    \n    :options => petsc_options(:ksplu),\n\n    :nls_trace =>true,\n    :nls_iter => 20,\n\n    :ν => 1.0e-5,  #channel = 0.0001472, \n    :ρ => 1.0, #kg/m3 density\n    :body_force => 0.0, #channel = 0.00337204\n    \n    :np_x => 2, #number of processors in X\n    :np_y => 2, #number of processors in Y\n    :np_z => 1, #number of processors in Z\n\n    :start_condition => :turbulent,\n    :restart => false,\n    :restart_file => \"DU89_0p85.csv\",\n    :TI =>0.01,\n    :Vbox => turbulence_box(),\n    :Re_filename => \"none\"\n      )\n\n    parameters = initialize_parameters(parameters)\n    \n    return parameters\nend\n\nparams = instantiate_parameters()\n\n\nExoFlow.main((params, backend))","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"The example above run in the REPL emulating a parallel run over 4 processors (you can see it by the options :np_x and :np_y).  The final results is a .pvd file, which can be open with ParaView. It is an index file of .pvtu files in Results/ folder. Each .pvtu file is the solution at a single time step and each .pvtu is and index file of .pvtu files. Each of this .vtu files is the portion of one processor of the solution of a spefic time step.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"info: numeric info\nChanging the backend to MPIBackend() allows it to run in MPI. ","category":"page"},{"location":"Cases/liddriven/#Lid-Driven-Cavity-Flow","page":"Lid Driven Cavity Flow","title":"Lid Driven Cavity Flow","text":"","category":"section"},{"location":"Cases/liddriven/","page":"Lid Driven Cavity Flow","title":"Lid Driven Cavity Flow","text":"(Image: Lid driven cavity flow) (Image: Lid driven cavity flow)","category":"page"},{"location":"Cases/liddriven/","page":"Lid Driven Cavity Flow","title":"Lid Driven Cavity Flow","text":"The Lid Driven Cavity flow is another standard case. It is a box of 1x1 dimension with the top side that can slide. The user can set different Reynolds.","category":"page"},{"location":"Cases/channel/#Channel","page":"Channel","title":"Channel","text":"","category":"section"},{"location":"Cases/channel/","page":"Channel","title":"Channel","text":"(Image: Re tau = 395)","category":"page"},{"location":"Cases/channel/","page":"Channel","title":"Channel","text":"The channel case can be run in 2D or 3D, with periodic boundaries or not. ","category":"page"},{"location":"Cases/airfoil/#Airfoil","page":"Airfoil","title":"Airfoil","text":"","category":"section"},{"location":"Cases/airfoil/","page":"Airfoil","title":"Airfoil","text":"(Image: SD7003 profile at Reynolds 60000)","category":"page"},{"location":"Cases/airfoil/","page":"Airfoil","title":"Airfoil","text":"It is one of the most complex and intersting case. The user has to create a proper mesh in gmsh setting the following physical boundaries:","category":"page"},{"location":"Cases/airfoil/","page":"Airfoil","title":"Airfoil","text":"inlet for the inlet\noutlet for the outlet\nairfoil for the airfoil walls\nlimits for the top and bottom boundaries","category":"page"},{"location":"Cases/airfoil/","page":"Airfoil","title":"Airfoil","text":"The velocity at the inlet is incresed from 0.0 arriving to the target value u_in at :t_endramp. This increase the numeric stability. If :t_endramp = :t0 the velocity at the inlet will be immediately :u_in. For numeric stability is better to keep u_in = 1.0, then fix the Reynolds and so the viscosity will be automatically computed as: ν = 1/Reynolds","category":"page"},{"location":"Cases/airfoil/","page":"Airfoil","title":"Airfoil","text":"The pressure is set 0.0 at the outlet section. The velocity on the limits is set equal to the one at inlet.","category":"page"},{"location":"Cases/cylinder/#Cylinder","page":"Cylinder","title":"Cylinder","text":"","category":"section"},{"location":"Cases/cylinder/","page":"Cylinder","title":"Cylinder","text":"(Image: Cylinder Vortex Shedding) (Image: Cylinder Vortex Shedding)","category":"page"},{"location":"Cases/cylinder/","page":"Cylinder","title":"Cylinder","text":"The Cylinder case can be used to see how meshes created in  gmsh are manged and to obtain the vortex shedding phenomena. The user has to create a proper mesh in gmsh setting the following physical boundaries:","category":"page"},{"location":"Cases/cylinder/","page":"Cylinder","title":"Cylinder","text":"inlet for the inlet\noutlet for the outlet\ncylinder for the cylinder walls\nlimits for the top and bottom boundaries","category":"page"},{"location":"numericsettings/#Numeric-Settings","page":"Numeric Settings","title":"Numeric Settings","text":"","category":"section"},{"location":"numericsettings/","page":"Numeric Settings","title":"Numeric Settings","text":"In each section, for each type of simulation, there are numerical hint and explation of how the boundaries and data are managed. ","category":"page"},{"location":"numericsettings/","page":"Numeric Settings","title":"Numeric Settings","text":"info: numeric info\nIn general, the body force should be zero. The only reasonable execption is for the periodic channel.","category":"page"},{"location":"numericsettings/","page":"Numeric Settings","title":"Numeric Settings","text":"info: numeric info\nIn general is better to set a reference speed u_in = 1 and the desired reynolds, the viscosity is then computed ad 1/Reynods. ","category":"page"},{"location":"numericsettings/","page":"Numeric Settings","title":"Numeric Settings","text":"warning: numeric waning\nIn MPI context the Additive Shwarz Method is shown to better scale than the Geometric Algebraic MultiGrid.","category":"page"},{"location":"numericsettings/","page":"Numeric Settings","title":"Numeric Settings","text":"info: numeric info\nIt is better to avoid using :julia as a solver because it uses a direct solver which is slow for large problems.","category":"page"},{"location":"ref/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"ref/#Conservation","page":"Reference","title":"Conservation","text":"","category":"section"},{"location":"ref/","page":"Reference","title":"Reference","text":"conservations_equations\nderivative_conservations_equations","category":"page"},{"location":"ref/#ExoFlow.conservations_equations","page":"Reference","title":"ExoFlow.conservations_equations","text":"conservations_equations(params::Dict{Symbol,Any})\n\nIt provides the Navier Stokes equations: continuity and momentum. The equations are different if the problem is linearized or not. Conitnuity: nablacdot u = 0\n\nMomentum: dfracpartial upartial t + u cdot nabla(u) + nabla (p) + nu Delta u= 0\n\nMomentum linearized: tildeu cdot nabla(u) + nabla (p) + nu Delta u = 0 where tildeu is the approximation of the convective velocity.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.derivative_conservations_equations","page":"Reference","title":"ExoFlow.derivative_conservations_equations","text":"derivative_conservations_equations(params::Dict{Symbol,Any})\n\nIt provides the spatial derivative of the Navier Stokes equations: continuity and momentum. The equations are different if the problem is linearized or not. Conitnuity: nablacdot du = 0\n\nMomentum: du cdot nabla(u) + u cdot nabla(du) + nabla (dp) + nu Delta du= 0\n\nMomentum linearized: tildeu cdot nabla(u) + nabla (p) + nu Delta u = 0\n\n\n\n\n\n","category":"function"},{"location":"ref/#Linear-Utilities","page":"Reference","title":"Linear Utilities","text":"","category":"section"},{"location":"ref/","page":"Reference","title":"Reference","text":"create_ũ_vector\nupdate_ũ\nupdate_ũ_vector!\nupdate_linear!","category":"page"},{"location":"ref/#ExoFlow.create_ũ_vector","page":"Reference","title":"ExoFlow.create_ũ_vector","text":"It instantiate an AbstactVector of 4 elements. Each element is an AbstractVector where the values of the velocity in each node is stored. The first element refers to the time step n, the second to the time step n-1, the third n-2 and the last n-3.\n\n\n\n\n\ncreateũvector(zfields::SequentialData)\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.update_ũ","page":"Reference","title":"ExoFlow.update_ũ","text":"It updates the convective velocity exitmation tildeu for the approximation of the non linear term: tildeu cdot nabla(u)\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.update_ũ_vector!","page":"Reference","title":"ExoFlow.update_ũ_vector!","text":"It updates the vector which stores the values of velocity at previous time steps.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.update_linear!","page":"Reference","title":"ExoFlow.update_linear!","text":"update_linear!(params::Dict{Symbol,Any})\n\nWrapper function for updating the ũ_vector containing the speed values at previous time steps  and ũ which is the approximation of ũ for the non linear terms, ũ⋅(∇u)\n\n\n\n\n\n","category":"function"},{"location":"ref/#SUPG","page":"Reference","title":"SUPG","text":"","category":"section"},{"location":"ref/","page":"Reference","title":"Reference","text":"h_param\nτ\nτb\nSUPG\nSUPG_lin","category":"page"},{"location":"ref/#ExoFlow.h_param","page":"Reference","title":"ExoFlow.h_param","text":"It computes the dimension of each cell of the computational domain. It uses the approximation: h = M^1D where h is the cell dimension, M the cell measure (surface for 2D, volume for 3D) and D is the Dimension (2 or 3)\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.τ","page":"Reference","title":"ExoFlow.τ","text":"Stabilization parameter for momentum equation for the SUPG.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.τb","page":"Reference","title":"ExoFlow.τb","text":"Stabilization parameter for continuity equation for the SUPG.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.SUPG","page":"Reference","title":"ExoFlow.SUPG","text":"SUPG non linear variational forumulation. It calles the conservations_equations and derivative_conservations_equations. It provides the equations set, the jacobian and the time jacobian.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.SUPG_lin","page":"Reference","title":"ExoFlow.SUPG_lin","text":"SUPG linear variational forumulation. It calles the conservations_equations and derivative_conservations_equations.  The terms are divided between mass matrix, stifness matrix and right hand side.\n\n\n\n\n\n","category":"function"},{"location":"ref/#VMS","page":"Reference","title":"VMS","text":"","category":"section"},{"location":"ref/","page":"Reference","title":"Reference","text":"G_params\nτm\nτc\nVMS\nVMS_lin","category":"page"},{"location":"ref/#ExoFlow.G_params","page":"Reference","title":"ExoFlow.G_params","text":"It provides the set of stabilization parameters for the VMS: G, GG, gg.  GG comes from the inverse of the jacobian of the map between reference and physical domain. It is exactly the cell dimension obtained by h_param for orthogonal carthesian grids.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.τm","page":"Reference","title":"ExoFlow.τm","text":"Stabilization parameter for momentum equation for the VMS.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.τc","page":"Reference","title":"ExoFlow.τc","text":"Stabilization parameter for contitnuity equation for the VMS.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.VMS","page":"Reference","title":"ExoFlow.VMS","text":"VMS non linear variational forumulation. It calles the conservations_equations and derivative_conservations_equations. It provides the equations set, the jacobian and the time jacobian.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.VMS_lin","page":"Reference","title":"ExoFlow.VMS_lin","text":"VMS linear variational forumulation. It calles the conservations_equations and derivative_conservations_equations.  The terms are divided between mass matrix, stifness matrix and right hand side.\n\n\n\n\n\n","category":"function"},{"location":"ref/#PETSc","page":"Reference","title":"PETSc","text":"","category":"section"},{"location":"ref/","page":"Reference","title":"Reference","text":"petsc_options","category":"page"},{"location":"ref/#ExoFlow.petsc_options","page":"Reference","title":"ExoFlow.petsc_options","text":"petsc_options(prec::Symbol)\n\nThe function has some solver options for PETSc solver. The user can set the argument of the function:\n\n- `:snesgamg` for using a Geometric Algebraic MultiGrid preconditioner and the Newton method for solving a non linear problem\n- `:kspgamg` for using a Geometric Algebraic MultiGrid preconditioner and the gmres method for solving a linear problem\n- `:sneslu` for using an Additive Shwarz Method as preconditioner and a LU factorization as sub-preconditioner and the Newton method for solving a non linear problem\n- `:ksplu` for using an Additive Shwarz Method as preconditioner and a LU factorization as sub-preconditioner and the gmres method for solving a linear problem\n\nThe gamg preconditioner is not compatible with the SequentialBackend()\n\n\n\n\n\n","category":"function"},{"location":"ref/#Common-procedures","page":"Reference","title":"Common procedures","text":"","category":"section"},{"location":"ref/","page":"Reference","title":"Reference","text":"hf_gen!\nadd_centre_tag!\nprintmodel\ncreation_fe_spaces\ncreation_op\ncreation_initial_conditions\ncreate_system_solver\ncreation_ode_parameters\ncreation_ode_solver\ncompute_solution_\ncompute_solution\niterate_solution\ncompute_solution_benchmark\nsolve_case\n\ninitialize_parameters","category":"page"},{"location":"ref/#ExoFlow.hf_gen!","page":"Reference","title":"ExoFlow.hf_gen!","text":"Add the body force. The user in the initial parameters specify the component in the x direction. Then the body force function is created 2D or 3D accordingly to the physics of the case.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.add_centre_tag!","page":"Reference","title":"ExoFlow.add_centre_tag!","text":"It creates the 'center' tag at the tag_coordinate (Point); if mesh extremely fine the tolrances have to be smaller (unlikely)\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.printmodel","page":"Reference","title":"ExoFlow.printmodel","text":"It prints the mesh of the model\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.creation_fe_spaces","page":"Reference","title":"ExoFlow.creation_fe_spaces","text":"It creates the finite elements spaces accordingly to the previously generated dirichelet tags\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.creation_op","page":"Reference","title":"ExoFlow.creation_op","text":"It creates the problem, here are implemented the VMS or SUPG system of equations\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.creation_initial_conditions","page":"Reference","title":"ExoFlow.creation_initial_conditions","text":"It creates the initial conditions. The default initial value for the pressure is 0.0 Pa. If the restart parameter is set to true, the initial conditons are set accordingly to the provided file\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.create_system_solver","page":"Reference","title":"ExoFlow.create_system_solver","text":"It creates the system solver which can be linear/non linear, julia/petsc using GAMG (only in petsc) or not.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.creation_ode_parameters","page":"Reference","title":"ExoFlow.creation_ode_parameters","text":"It creates parameters that are used for the ODE solution.  In specific, if different theta are specified for velocity and pressure, are created vectors which relate each degree of freedom to the desired theta .\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.creation_ode_solver","page":"Reference","title":"ExoFlow.creation_ode_solver","text":"It creates the ODE solver using the parameters defined by the user\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.compute_solution_","page":"Reference","title":"ExoFlow.compute_solution_","text":"It computes the solution, different for benchmark case or for getting the results.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.compute_solution","page":"Reference","title":"ExoFlow.compute_solution","text":"It computes the solution, in the serial or parallel case.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.iterate_solution","page":"Reference","title":"ExoFlow.iterate_solution","text":"It iterates the solution, here is where the solution is actually computed. At the end of each step the forces values are extracted if required\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.compute_solution_benchmark","page":"Reference","title":"ExoFlow.compute_solution_benchmark","text":"For benchmarking. It removes the compilation time for the first iteration and it does not print any results.\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.solve_case","page":"Reference","title":"ExoFlow.solve_case","text":"Wrapper function that solve the case. Here the Non Linear Solver - nls is defined, in case the default julia sovler is adopted or PETSc\n\n\n\n\n\n","category":"function"},{"location":"ref/#ExoFlow.initialize_parameters","page":"Reference","title":"ExoFlow.initialize_parameters","text":"initialize_parameters(params::Dict{Symbol, Any})\n\nIt makes some checks on the validity of the input provided by the user and add some additional information.\n\n\n\n\n\n","category":"function"},{"location":"#ExoFlow.jl","page":"Introduction","title":"ExoFlow.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Documentation of ExoFlow.jl for solving incompressible Navier-Stokes using stabilized Finite Element Method, in specific Streamline-Upwind Petrov-Galerkin (SUPG) and Variational Multiscale Method (VMS)","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"At the moment it is not a registered Julia package. For installing, from the REPL just press ].","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(@1.8) pkg> add https://github.com/carlodev/ExoFlow.jl ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"It uses some custom forks of registered julia packages which need to be installed manually. For installing, from the REPL just press ].","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(@1.8) pkg> add https://github.com/carlodev/Gridap.jl#Thetamethod_multifield\n(@1.8) pkg> add https://github.com/carlodev/GridapDistributed.jl#Thetamethod_multifield\n(@1.8) pkg> add https://github.com/carlodev/SyntheticEddyMethod.jl#master\n(@1.8) pkg> add https://github.com/carlodev/GridapPETSc.jl#master","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"For a complete and smooth experience is suggested to install also the free software ParaView which allows to graphically visualize the results.","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"#Package-General-Features","page":"Introduction","title":"Package General Features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Implementation of SUPG and VMS formulation for same order elements for velocity and pressure\nUse of different ode method: θ-method with different values of θ for velocity and pressure, and α-mehtod\nUsing custom Meshes created with gmsh. For airfoils the package AirfoilGmsh.jl has been developed for speeding up the process\nSolve 2D and 3D cases\nPossibility of chosing the backend thanks to PartitionedArrays.jl. It can be run in the REPL or in MPI","category":"page"},{"location":"#Known-Issues","page":"Introduction","title":"Known Issues","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The gamg preoconditioner is not compatible with SequentialBackend()\nThe final results is a .vtu file, which can be open with ParaView. It is an index file of .pvtu files in Results/ folder. Each .pvtu file is the solution at a single time step and each .pvtu is and index file of .pvtu files. Each of this .vtu files is the portion of one processor of the solution of a spefic time step. The .pvtu files have a wrong indexation, the user has to remove the Results/.","category":"page"},{"location":"#Acknowledgement","page":"Introduction","title":"Acknowledgement","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The Variational Multiscale formualation is based on the paper: Bazilevs Y. et al., Variational multiscale residual-based turbulence modeling forlarge eddy simulation of incompressible flows , 10.1016/j.cma.2007.07.016","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The Streamline Upwind Petrov Galerkin is based on the PhD thesis:  Tamas Banyai,  Development of Stabilized Finite ElementMethod for Numerical Simulation ofTurbulent Incompressible Single andEulerian-Eulerian Two-Phase Flows ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The linearization of the equations is based on the PhD thesis:  Bart Janssen, Numerical modeling and experimentalinvestigation of fine particle coagulationand dispersion in dilute flows ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The package is based on the use of Gridap: Santiago Badia and Francesc Verdugo,  Gridap: An extensible Finite Element toolbox in Julia, 10.21105/joss.02520","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The solver used in MPI is PETSc: Satish Balay et al., https://petsc.org/","category":"page"},{"location":"Cases/taylorgreen/#Taylor-Green-Vortex","page":"Taylor Green","title":"Taylor Green Vortex","text":"","category":"section"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"(Image: Taylor Green) (Image: Taylor Green)","category":"page"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"It solves the 2D Taylor Green Vortex case. It is the only case where analtical solution are available:","category":"page"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"analytical_solution","category":"page"},{"location":"Cases/taylorgreen/#ExoFlow.analytical_solution","page":"Taylor Green","title":"ExoFlow.analytical_solution","text":"analytical_solution(diameter::Int64, Vs::Float64, Ua::Float64, Va::Float64, ν::Float64) \n\nIt provides the anlytical solution for the Taylor Green Vortex case.\n\n\n\n\n\n","category":"function"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"u_x= U_a-V_s cos bigg (fracpiD(x-U_a t)bigg ) sin bigg (fracpiD(y-V_a t)bigg ) e^-frac2 v pi^2D^2 t","category":"page"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"u_y= V_a+V_s sin bigg (fracpiD(x-U_a t)bigg ) cos bigg (fracpiD(y-V_a t)bigg ) e^-frac2 nu pi^2D^2 t","category":"page"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"p=-fracV_s^24bigg ((cos(2 fracpiD(x-U_a t) )+cos (2 fracpiD(y-V_a t))bigg ) e^-frac4 v pi^2D^2 t","category":"page"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"omega=frac2 V_s piD cos bigg (fracpiD(x-U_a t)bigg ) cos bigg (fracpiD(y-V_a t)bigg ) e^-frac4 nu pi^2D^2 t","category":"page"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"Because of an analytical solution it is used as a benchamark case for verifing mesh convergence, CFL stability and processor scalability. The domain is a squqare of 2Dx2D with periodic boundaries over the 4 sides. The initial solution is retrived from the analytical solution. The pressure is fixed in the centre of the domain equal to the analytical solution.  The parameters set by the user are overwritten by the following standard values:","category":"page"},{"location":"Cases/taylorgreen/","page":"Taylor Green","title":"Taylor Green","text":"Vs = 1ms swirling velocity\nUa = 02ms translational velocity in x direction\nVa = 03ms translational  velocity in y direction\nD = 05m vortex dimensions\nnu = 0001 m^2s","category":"page"}]
}
