function main(main_args::Tuple)
    params, backend = main_args
    with_backend(backend, params[:procs_partition]) do parts

    if backend == MPIBackend()
        comm = MPI.COMM_WORLD
        #To avoid multiple printing of the same line in parallel
        if MPI.Comm_rank(comm) != 0
            redirect_stderr(devnull)
            redirect_stdout(devnull)
        end    
        params = merge!(params, Dict(:comm => comm))
    end
 
    params = merge!(params, Dict(:parts => parts, :backend => backend))
    

    # params = instantiate_parameters()
    params = merge!(params, Dict(:parts => parts))
    println(params)
    #Unpack the dictionary
    
    @unpack case, solver, options = params
    if case == "Channel"
        if solver == :petsc
            GridapPETSc.with(args=split(options)) do
                
                run_channel(params)
            end
        elseif solver == :julia
            run_channel(params)
        end

    elseif case == "LidDriven"
        if  solver == :petsc
            GridapPETSc.with(args=split(options)) do
                run_liddriven(params)
            end
        elseif solver == :julia
            run_liddriven(params)
        end


    elseif case == "Cylinder"
        if  solver == :petsc
            GridapPETSc.with(args=split(options)) do
                run_cylinder(params)
            end
        elseif solver == :julia
            run_cylinder(params)
        end

    elseif case == "Airfoil"
        if  solver == :petsc
            GridapPETSc.with(args=split(options)) do
                run_airfoil(params)
            end
        elseif solver == :julia
            run_airfoil(params)
    end

    elseif case == "TaylorGreen"
        if  solver == :petsc
            GridapPETSc.with(args=split(options)) do
                run_taylorgreen(params)
            end
        elseif solver == :julia
            run_taylorgreen(params)
    end

elseif case == "Plate"
    if  solver == :petsc
        GridapPETSc.with(args=split(options)) do
            run_plate(params)
        end
    elseif solver == :julia
        run_plate(params)
end

elseif case == "Jet"
        #run_jet(parts; info)


    else
        error("Case: $(case) not defined")
    end
    
    #end
end


end 