function extract_matrix_vector(params)
    if params[:ode_method] == :AlphaMethod
        error("Debug mode, extracting_matrix_vector only with ThetaMethod")       
    end
    @unpack xh0, op, ode_solver, dt, t0 = params
    
    uh0 = xh0

 
    u00 = get_free_dof_values(uh0)


 
    
    uf = copy(get_free_dof_values(uh0))
    odeop = get_algebraic_operator(op)
 

    ode_cache = allocate_cache(odeop)
    vθ = similar(u00)
  
    nl_cache = nothing

    ode_solver.θ == 0.0 ? dtθ = dt : dtθ = dt*ode_solver.θ
   
    tθ = t0+dtθ
  
    ode_cache = update_cache!(ode_cache, odeop, tθ)
    println("extracting matrix and vector")
    nlop = Gridap.ODEs.ODETools.ThetaMethodNonlinearOperator(odeop,tθ,dtθ,u00,ode_cache,vθ)
    

    nl_cache = solve!(uf, ode_solver.nls, nlop, nl_cache)
    
    Anl = nl_cache.f0
    bnl = nl_cache.j0
    
    comm = MPI.COMM_WORLD
    for i = 0:1:MPI.Comm_size(comm)-1
      if (MPI.Comm_rank(comm)) == i
      io = open("output_vector_matrix_$i.txt", "w")
      write(io, "Anl = $(Anl.values) \n")
      write(io, "bnl = $(bnl.values) \n")   
      close(io)
    end
    println("Partitioned Matrices printed in files")
    end
    return bnl.values, Anl.values
  end

  function print2file_stab_params(params, stab_params)
 
    comm = MPI.COMM_WORLD
   
    if params[:debug_mode]
      if params[:method] == :VMS
       stab_params_names = ["G", "GG", "gg"]
      
      elseif params[:method] == :SUPG
        stab_params_names = ["h"]        
      end

      for i = 0:1:MPI.Comm_size(comm)-1
        
        if (MPI.Comm_rank(comm)) == i
        filename  = Dates.format(now(), "HMSs") * "stabilization_parameter_$(i).txt"
        println(filename)
        io = open(filename , "w")
        println("file opening")
        for k = 1:1:length(stab_params_names)
          println("file writing")
          println(io, "$(stab_params_names[k]) = $(stab_params[k]) \n")
        end
   
        close(io)
      
      end

      println("Stabilization parameters printed in files")
      end
      MPI.Barrier(comm)
    end
  
  end

  "Use a low level custom assembler for petsc"
  function PETSC_TransientFEOperator(res, jac, jac_t, X, Y)
    Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
    Tv = Vector{PetscScalar}
    assem_t=SparseMatrixAssembler(Tm, Tv, X(nothing), Y)
    op = Gridap.ODEs.TransientFETools.TransientFEOperatorFromWeakForm{Nonlinear}(res,(jac,jac_t),assem_t,(X,∂t(X)),Y,1)
   
    return op
  end