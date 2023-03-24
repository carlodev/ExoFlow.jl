using ExoFlow
using PartitionedArrays


include("Turbulence_Settings.jl")


#channel body force = 0.00337204, ν=0.0001472, u_in =1.0
#lid driven u_in = 0

function driver_test(case, method, odemethod; linear = true, D = 2 )
if case == "Airfoil"
mesh_gen = true
else
    mesh_gen = false
end

if linear
    petsc_opt = :ksplu
else
    petsc_opt = :sneslu
end

backend = SequentialBackend() #or MPIBackend() SequentialBackend()
function instantiate_parameters()
    parameters = Dict(
    :N => 20,
    :D => D, #Dimension
    :order => 1,
    
    :t0 => 0.0,
    :dt => 0.01,
    :tF => 0.05,
    :t_endramp => 0.0,

    :case => case,
    :solver => :petsc,
    :method => method,
    :ode_method => odemethod,
    :θ=>[0.5, 1.0],
    :ρ∞ => 1.0,
    :Re => 1_000,
    :c => 1, #chord lenght [m], used for naca (=1), cylinder (=0.1), liddriven = (1), 0.5
    :u_in => 1.0,  # =1.0 for lid driven 
    :periodic => false,

    :printmodel => false,
    
    :mesh_gen => mesh_gen,

    :linear => linear,
    :steady => false,

    :debug_mode => false,
    :benchmark_mode => true,
    :printinitial => false,

    :mesh_file => "DU89_precompile_2D.msh",
    :Cᵢ => [4, 36],    
    :options => petsc_options(petsc_opt),

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
      )

    parameters = initialize_parameters(parameters)
    
    return parameters
end

params = instantiate_parameters()


ExoFlow.main((params, backend))

return true
end