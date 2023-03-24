include("MeshChannel.jl")
include("BoundaryConditions.jl")




function run_channel(params)

  
    #model_serial = mesh_channel(nothing, D, N, false, periodic)
    model = mesh_channel(params)
    println("Model created")
    printmodel(params, model)
    
    u0 = initial_velocity(params)
    println("u0 created")
    merge!(params, Dict(:u0 => u0, :model => model))

   hf_gen!(params)

    u_diri_tags, u_diri_values, p_diri_tags, p_diri_values = bc_channel(params)
    println("tags created")
    
    
    V, U, P, Q, Y, X = creation_fe_spaces(params, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)
    println("Spaces created")


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


end