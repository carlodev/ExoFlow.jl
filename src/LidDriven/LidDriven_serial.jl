include("../ExoFlow.jl")
include("BoundaryConditions.jl")

function petsc_gamg_options()
    "-snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 1.0 -snes_rtol 1.0e-14 -snes_atol 0.0 -snes_monitor -pc_type asm -sub_pc_type lu -ksp_type gmres -ksp_gmres_restart 30  -snes_converged_reason -ksp_converged_reason -ksp_error_if_not_converged true"
end

    function stretching_y_function(x)
        gamma1 = 2.5
        S = 0.5815356159649889 #for rescaling the function over the domain -0.5 -> 0.5
        -tanh.(gamma1 .* (x)) ./ tanh.(gamma1) .* S
    end

    
    function stretching(x::Point)
        m = zeros(length(x))
        m[1] = stretching_y_function(x[1])


        m[2] = stretching_y_function(x[2])
        Point(m)
    end



    parameters = Dict(
        :N => 10,
        :D => 2, #Dimension
        :order => 1, 
        :t0 =>0.0,
        :dt => 0.5, #0.01
        :tF => 500,
        :t_endramp => 0.0,
    
        :case => "LidDriven",
        :solver => :petsc,
        :method => :VMS,
        :ode_method => :ThetaMethod, 
        :ρ∞ => 0.8, #0.01
        :θ=>[0.5,1.0],
        :Re => 10_000,
        :c => 1.0, #chord lenght [m], used for naca (=1), cylinder (=0.1), liddriven = 1
        :u_in => 1.0,  # =1.0 for lid driven 
        :linear =>false,
        :periodic => false,
        :printmodel => true,

        :mesh_gen => false,

        :linear => false,
        :steady => false,
    
        :debug_mode => false,
        :benchmark_mode =>false,
        
        :mesh_file => "NACA0012_2D_improved.msh",
        :Cᵢ => [10, 1], 
        :options => petsc_gamg_options(), 
        :nls_trace =>true,
        
        :ν => 4.61e-5, #channel = 0.0001472, 
        :ρ => 1.0, #kg/m3 density
        :body_force => 0.0, #channel = 0.00337204
    
        :np_x => 8, #number of processors in X
        :np_y => 8, #number of processors in Y
        :np_z => 2, #number of processors in Z
    
        :restart => false,
        :restart_file => "Airfoil/DU89_319.csv"
    ) 
    params = initialize_parameters(parameters)

    if params[:mesh_gen]
        mesh_file_path = joinpath(@__DIR__, "../../models", params[:mesh_file])

        model = GmshDiscreteModel( mesh_file_path)        
    else
        L = 0.5
        domain = (-L, L, -L, L)
        partition = (params[:N], params[:N])
        model = CartesianDiscreteModel(domain, partition, map=stretching)
        
    end

    printmodel(params, model)
    hf_gen!(params)

    u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u0 = bc_liddriven(model, params) #u_top in bc_liddriven taken as initial velocity
    merge!(params, Dict(:u0 => u0, :model => model))


    
    
     V, U, P, Q, Y, X = creation_fe_spaces(params, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)
    println("spaces created")

    degree = 4*params[:order]
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    new_dict = Dict(:U => U,
                    :P => P,
                    :X => X,
                    :Y => Y,
                    :Ω => Ω,
                    :dΩ => dΩ,
                    :degree => degree,
                    :force_params => nothing)
    merge!(params, new_dict)

    
    solve_case(params)


   

    

  

