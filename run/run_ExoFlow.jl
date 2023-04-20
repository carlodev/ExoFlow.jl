using ExoFlow
using PartitionedArrays
include("Turbulence_Settings.jl")



#channel body force = 0.00337204, ν=0.0001472, u_in =1.0
#lid driven u_in = 0


backend = MPIBackend() #or MPIBackend() SequentialBackend()
function instantiate_parameters()
    parameters = Dict(
    :N => 10,
    :D => 2, #Dimension
    :order => 1,
    
    :t0 => 0.0,
    :dt => 0.0015,
    :tF => 10.0,
    :t_endramp => 0.0,

    :case => "Plate",
    :solver => :petsc,
    :method => :VMS,
    :ode_method => :ThetaMethod,
    :θ=>[0.5, 1.0],
    :ρ∞ => 1.0,
    :Re => 1_000,
    :c => 1, #chord lenght [m], used for naca (=1), cylinder (=0.1), liddriven = (1), 0.5
    :u_in => 1.0,  # =1.0 for lid driven 
    :periodic => true,

    :printmodel => true,
    
    :mesh_gen => false,

    :linear => true,
    :steady => false,

    :debug_mode => false,
    :benchmark_mode => false,
    :printinitial => true,
    :print_last =>false, 
  
    :mesh_file => "sd7003s_3D.msh",
    :Cᵢ => [4, 36],    
    :options => petsc_options(:ksplu),

    :nls_trace =>true,
    :nls_iter => 20,

    :ν => 1.0e-5,  #channel = 0.0001472, 
    :ρ => 1.0, #kg/m3 density
    :body_force =>0.0, #channel = 0.00337204
    
    :np_x => 2, #number of processors in X
    :np_y => 2, #number of processors in Y
    :np_z => 1, #number of processors in Z

    :start_condition => :turbulent,
    :restart => false,
    :restart_file => "src/Airfoil/Sd1.csv",
    :TI =>0.01,
    :Vbox => turbulence_box(),
    :Re_filename => "none"
      )

    parameters = initialize_parameters(parameters)
    
    return parameters
end

params = instantiate_parameters()


ExoFlow.main((params, backend))

#mpiexecjl --project=../ -n 4 julia -O3 --check-bounds=no -L run_ExoFlow.jl