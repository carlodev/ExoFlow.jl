include("BoundaryConditions.jl")



function run_liddriven(params)
    

    function stretching_y_function(x)
        gamma1 = 2.5
        S = 0.5815356159649889 #for rescaling the function over the domain -0.5 -> 0.5
        -tanh.(gamma1 .* (x)) ./ tanh.(gamma1) .* S
    end

    
    function stretching(x::Point)
        m = zeros(length(x))
        m[1] = stretching_y_function(x[1])
        m[2] = stretching_y_function(x[2])
        if length(x)>2
            m[3] = stretching_y_function(x[3])
        end
        Point(m)
    end

    if params[:mesh_gen]
        mesh_file_path = joinpath(@__DIR__, "../../models", params[:mesh_file])

        model = GmshDiscreteModel(params[:parts], mesh_file_path)

    else
        L = 0.5
        if params[:D] == 2
            domain = (-L, L, -L, L)
            partition = (params[:N], params[:N])
        elseif params[:D] == 3
            println("Model 3 a")
            domain = (-L, L, -L, L, -L, L)
            partition = (params[:N], params[:N], params[:N])
        end

    
        model = CartesianDiscreteModel(params[:parts], domain, partition, map=stretching)
        
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


   

    

  

end