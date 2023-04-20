include("BoundaryConditions.jl")
include("Blasius.jl")
"""
Mesh and dimension from:
Thomas J R Hughes and Victor M Calo and Guglielmo Scovazzi, VARIATIONAL AND MULTISCALE METHODS IN TURBULENCE, 2005

ν = 2e-6
δ0 = 1.59e-3
xin = 0.005 
xtr = 1
"""


function run_plate(params)

  params[:ν] = 2e-5 #m2/s 
  params[:u_in] = 1.0
  δ0 = 1.59e-3
  xin = 0.005 

  merge!(params, Dict(:δ0 =>δ0, :xin=>xin))
  Lx = 620*δ0
  Ly = 40*δ0
  Lz = 30*δ0

  Nx = 512
  Ny = 48
  Nz = 48


    # dim_not_supported(params[:D], params[:case])
    function stretching(x::Point)
      D = length(x)
      m = zeros(D)
  if D == 2
       
        m[1] = x[1]
        m[2] = x[2].^3 / Ly^2
  elseif D == 3
        m[1] = x[1]
        m[2] = x[2].^3 / Ly^2
        m[3] = x[3]
      end
        Point(m)
  end
      

 
    @unpack N,D = params
    if D == 2
    
      domain = (0, Lx, 0, Ly)
      partition = (Nx, Ny)
      periodic_tuple = (false, false)
    elseif D == 3
      domain = (0, Lx, 0, Ly, -Lz/2, Lz/2)
      partition = (Nx, Ny, Nz)
      periodic_tuple = (false, false,true)
    end
    
    model = CartesianDiscreteModel(params[:parts], domain, partition,map=stretching, isperiodic = periodic_tuple)


  

 
    
    printmodel(params, model)
    hf_gen!(params)

    u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u0 = bc_plate(params) #u_top in bc_liddriven taken as initial velocity
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

