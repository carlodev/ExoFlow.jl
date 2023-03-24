include("../ExoFlow.jl")
include("SpaceConditions.jl")
include("AnalyticalSolution.jl")
# function petsc_options()
#   "-snes_type newtonls -snes_linesearch_type basic -snes_linesearch_damping 1.2 -snes_rtol 0.0 -snes_stol 1.0e-10 -snes_atol 0.0 -snes_monitor -pc_type asm -sub_pc_type lu   -ksp_type pgmres -ksp_gmres_restart 30  -snes_converged_reason -ksp_converged_reason -ksp_error_if_not_converged true"
# end

params = Dict(
  :N => 32,
  :D => 2, #Dimension
  :order => 1,
  
  :t0 => 0.0,
  :dt => 0.01,
  :tF => 0.5,
  :t_endramp => 1.0,


  :case => "TaylorGreen",
  :solver => :petsc,
  :method => :VMS,
  :ode_method => :ThetaMethod,
  :θ=>0.5,
  :ρ∞ => 0.8,
  :Re => 100_000,
  :c => 1, #chord lenght [m], used for naca (=1), cylinder (=0.1), liddriven = (1), 0.5
  :u_in => 1.0,  # =1.0 for lid driven 
  :periodic => false,

  :printmodel => false,
  
  :mesh_gen => false,

  :linear => false,
  :steady => false,

  :debug_mode => false,

  :mesh_file => "NACA0012_2D_improved.msh",
  :Cᵢ => [4, 36],    
  :options => petsc_options(),
  :nls_trace =>true,
  :nls_iter => 20,

  :ν => 1.0e-5,  #channel = 0.0001472, 
  :ρ => 1.0, #kg/m3 density
  :body_force => 0.0, #channel = 0.00337204
  
  :np_x => 2, #number of processors in X
  :np_y => 2, #number of processors in Y
  :np_z => 1, #number of processors in Z

  :restart => false,
  :restart_file => "Du89_2p1.csv",
  :TI =>0.01,

  )
  
  params= initialize_parameters(params)
  diameter = 0.5 #0.5 [m] vortex dimension

  Vs = 1 #1[m/s]swirling speed
  Ua = 0.3 #0.3 [m/s]convective velocity in x
  Va = 0.2 #0.2 [m/s]convective velocity in y
  params[:ν] = 0.001 #0.001 m2/s 
  #Domain and mesh definition



    #Domain and mesh definition
    domain = (-diameter, diameter, -diameter, diameter)
    partition = (params[:N], params[:N])
    model = CartesianDiscreteModel(domain, partition; isperiodic=(true, true))
  
    hf_gen!(params)
  
    velocity, pa, ωa = analytical_solution(diameter, Vs, Ua, Va, params[:ν])
    merge!(params, Dict(:u0 => velocity, :model => model))
    V, Q, U, P, Y, X, model = CreateTGSpaces(model, params, pa) #We update model with the new label of the center point
    println("spaces created")
  
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
    
    xh0 = creation_initial_conditions(params)
    println("initial conditions created")
    merge!(params, Dict(:xh0 => xh0)) 
    ph0 =  interpolate_everywhere(pa(0), P(0))
    
    ph0_free_dofs = get_free_dof_values(ph0)
    
    pfun = FEFunction(Q, ph0_free_dofs )

    writevtk(Ω, "TG.vtu", cellfields=["ph0"=>pfun])

    vec_ph0 = [ph0_free_dofs, ph0_free_dofs, ph0_free_dofs, ph0_free_dofs]
  
    p_dofs = length(get_free_dof_values(ph0))  
  
    
    oscill_coeff = [-5, 9, -5, 1]
    function compute_stabilization(vec_ph0, ph0, oscill_coeff)
      Δph0 = ph0 - vec_ph0[1]  
      Δpc = 0.5*(2*Δph0 + oscill_coeff[1] * vec_ph0[1] + oscill_coeff[2] * vec_ph0[2]+ oscill_coeff[3] * vec_ph0[3] + oscill_coeff[4] * vec_ph0[4])
      return  ph0 - Δpc
    end
    
    function  update_sol_vec!(vec_ph0, ph0)
      vec_ph0[2:end] =  vec_ph0[1:end-1]
      vec_ph0[1] = ph0
    end
  
     
    cache = GridapPETSc.with(args =split(params[:options])) do
    op = creation_op(params)
    println("Op created")
    merge!(params, Dict(:op => op))



    if params[:solver] == :julia
      if params[:linear]
          sys_solver = LUSolver()
      else
        println("Julia non Linear Solver")
          sys_solver = NLSolver(show_trace=params[:nls_trace], method=:newton, iterations=params[:nls_iter])
      end
  
  elseif params[:solver] == :petsc
      if params[:linear]
          sys_solver = PETScLinearSolver()
          println("Using PETSc linear solver")
      elseif params[:petsc_snes] 
          sys_solver  =  PETScNonlinearSolver()
          println("Using PETSc non linear solver")
      elseif params[:petsc_snes] == false
          ls = PETScLinearSolver()
          println("Using PETSc linear solver")
          sys_solver = NLSolver(ls, show_trace=params[:nls_trace], method=:newton, iterations=params[:nls_iter])
          println("Using Gridap non linear solver")
      end
     
  end

  merge!(params, Dict(:sys_solver => sys_solver))
  println("system solver created")
  
  



  sol_t, ode_solver = creation_ode_solver(params)
  println("sol_t created")
  merge!(params, Dict(:sol_t => sol_t, :ode_solver => ode_solver))
  


  current, state = Base.iterate(sol_t)


  uf, tf = current

  while tf < params[:tF]
  
  current, state = Base.iterate(sol_t, state)
  uf, tf = current
  Uh, odesolstate = state
  uf_s,u0_s,t0,cache = odesolstate
  ph_tn = uf_s[end-p_dofs+1:end]
  
  uh_tn = uf[1]
  p_stab = compute_stabilization(vec_ph0, ph_tn, oscill_coeff)
  update_sol_vec!(vec_ph0, ph_tn)
  uf_s[end-p_dofs+1:end] .=p_stab
  uf_s[end-p_dofs+1:end] .=p_stab
  println(length(p_stab))

 odesolstate =   (uf_s,u0_s,t0,cache)
  state =  (Uh, odesolstate)
  p_analytic = pa(tf)
  ph_f = FEFunction(Q, p_stab)
  writevtk(Ω, "Results/LidDriven_$tf.vtu", cellfields=["uh" =>uf[1], "ph" =>uf[2], "ph_stab" =>ph_f, "p_analytic"=>p_analytic])
  println(tf)
  end

end

  
    #This is a "special case", the only one with analytical solutions
    e_u = 10
    e_p = 10
    iter = 0
    createpvd("$(params[:case])_$(params[:D])d") do pvd
      
      @time for (xh_tn, tn) in params[:sol_t]
    
        uh_tn = xh_tn[1]
        ph_tn = xh_tn[2]
        ωh_tn = ∇ × uh_tn
  
        Δu = Δ(uh_tn)
        p_analytic = pa(tn)
        u_analytic = velocity(tn)
        w_analytic = ωa(tn)
     
        if isapprox(tn, params[:tF]; atol=0.2 * params[:dt])
          e_u
          e_p
          eu = velocity(params[:tF]) - uh_tn
          ep = pa(params[:tF]) - ph_tn
          e_u = sqrt(sum(∫(eu ⋅ eu) * params[:dΩ]))
          e_p = sqrt(sum(∫(ep * ep) * params[:dΩ]))
  
        end #end if
        iter = iter + 1 
        println("iter = $iter")
        if params[:linear]
          println("Linear")
          if iter == 1
              global ũ_vector = [uh_tn, uh_tn, uh_tn, uh_tn] 
              global ũ = uh_tn
          else
          update_ũ_vector!(ũ_vector, uh_tn)
          global ũ = update_ũ(ũ_vector, params[:ũ_coeff])
          end
      end
      

        pvd[tn] = createvtk(params[:Ω], "Results/$(params[:case])_$tn" * ".vtu", cellfields=["uh" => uh_tn, "ph" => ph_tn, "ph_stab"=>ph_tn_stab, "wh" => ωh_tn, "p_analytic"=>p_analytic, "u_analytic"=>u_analytic,  "w_analytic"=>w_analytic])
      end #end for
  
    end #end do


e_u
e_p


err = [e_u, e_p] #error on velocity and on pressure
println("err_u =$(e_u)")
println("err_p =$(e_p)")
jldsave("$(params[:case])_error_$(params[:N]).jld2"; err)


