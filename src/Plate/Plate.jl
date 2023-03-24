include("BoundaryConditions.jl")



function run_plate(params)
    dim_not_supported(params[:D], params[:case])
    function stretching(x::Point)
        m = zeros(length(x))
        m[1] = x[1]
        m[2] = x[2].^2
        Point(m)
      end
      

    params[:ν] = 1e-5 #m2/s 
    params[:u_in] = 10.0

    Lx = 3
    Ly = 1
    Nx = 250
    Ny = 400
    

    domain = (0, Lx, 0, Ly)
    partition = (Nx, Ny)
    model = CartesianDiscreteModel(params[:parts], domain, partition,map=stretching)
    
    printmodel(params, model)
    hf_gen!(params)

    u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u0 = bc_plate(model, params[:u_in]) #u_top in bc_liddriven taken as initial velocity
    merge!(params, Dict(:u0 => u0, :model => model))


    
    
    U,P,Y,X = creation_fe_spaces(params, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)
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