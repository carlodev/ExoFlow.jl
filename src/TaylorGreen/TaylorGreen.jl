include("SpaceConditions.jl")
include("AnalyticalSolution.jl")


function run_taylorgreen(params)

  diameter = 0.5 #0.5 [m] vortex dimension

  Vs = 1 #1[m/s]swirling speed
  Ua = 0.3 #0.3 [m/s]convective velocity in x
  Va = 0.2 #0.2 [m/s]convective velocity in y
  params[:ν] = 0.001 #0.001 m2/s 

  #Domain and mesh definition
  domain = (-diameter, diameter, -diameter, diameter)
  partition = (params[:N], params[:N])
  model = CartesianDiscreteModel(params[:parts], domain, partition; isperiodic=(true, true))

  hf_gen!(params)

  velocity, pa, ωa = analytical_solution(diameter, Vs, Ua, Va, params[:ν])
  merge!(params, Dict(:u0 => velocity, :model => model))
  V, Q, U, P, Y, X, model = CreateTGSpaces(model, params, pa) #We update model with the new label of the center point
  @info "spaces created"


  printmodel(params, model)




  degree = 4 * params[:order]
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)

  new_dict = Dict(:U => U,
    :P => P,
    :X => X,
    :Y => Y,
    :Ω => Ω,
    :dΩ => dΩ,
    :degree => degree,
    :force_params => nothing,
    :p0 => pa)

  merge!(params, new_dict)



  sys_solver = create_system_solver(params)
  @info "system solver created"
  merge!(params, Dict(:sys_solver => sys_solver))

  xh0 = creation_initial_conditions(params)
  @info "initial conditions created"
  merge!(params, Dict(:xh0 => xh0)) 
  
  op = creation_op(params)
  @info "Op created"
  merge!(params, Dict(:op => op))

  sol_t, ode_solver = creation_ode_solver(params)
  @info "sol_t created"
  merge!(params, Dict(:sol_t => sol_t, :ode_solver => ode_solver))

  if params[:benchmark_mode]
    compute_solution_benchmark(params)
  
  else
    #This is a "special case", the only one with analytical solutions
    e_u = 10
    e_p = 10
    iter = 0
   @time createpvd(params[:parts], "$(params[:case])_$(params[:D])d") do pvd
      
       for (xh_tn, tn) in params[:sol_t]
    
        uh_tn = xh_tn[1]
        ph_tn = xh_tn[2]
        ωh_tn = ∇ × uh_tn
  
        Δu = Δ(uh_tn)
        p_analytic = pa(tn)
        u_analytic = velocity(tn)
        w_analytic = ωa(tn)
  
        # if isapprox(tn, params[:tF]; atol=0.2 * params[:dt])
        #   e_u
        #   e_p
        #   eu = velocity(params[:tF]) - uh_tn
        #   ep = pa(params[:tF]) - ph_tn
        #   e_u = sqrt(sum(∫(eu ⋅ eu) * params[:dΩ]))
        #   e_p = sqrt(sum(∫(ep * ep) * params[:dΩ]))
  
        # end #end if
        iter = iter + 1 
        println("iter = $iter")
        
        pvd[tn] = createvtk(params[:Ω], "Results/$(params[:case])_$tn" * ".vtu", cellfields=["uh" => uh_tn, "ph" => ph_tn, "wh" => ωh_tn, "p_analytic"=>p_analytic, "u_analytic"=>u_analytic,  "u0"=>params[:ũ]])
        
        
        update_linear!(params, uh_tn)

      
      
      end #end for
  
    end #end do
  end
  #   comm = MPI.COMM_WORLD
  
  #   if ((MPI.Comm_rank(comm)) == 0)
  #     err = [e_u, e_p] #error on velocity and on pressure
  #     println("err_u =$(e_u)")
  #     println("err_p =$(e_p)")
  
  #     jldsave("$(params[:case])_error_$(params[:N]).jld2"; err)
  #   end
  
  end