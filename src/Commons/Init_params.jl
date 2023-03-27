
"""
    initialize_parameters(params::Dict{Symbol, Any})

It makes some checks on the validity of the input provided by the user and add some additional information.
"""
function initialize_parameters(params::Dict{Symbol, Any})
    #processor partition, different if the mesh is internally generated or imported from a .msh file
     if params[:mesh_gen]
         if params[:D] == 2
             procs_partition = params[:np_x]* params[:np_y]
         elseif params[:D] == 3
             procs_partition = params[:np_x]*params[:np_y]* params[:np_z]
         else
             dim_error(params[:D])
         end

     else
        if params[:D] == 2
            procs_partition = (params[:np_x], params[:np_y])
        elseif params[:D] == 3
            procs_partition = (params[:np_x], params[:np_y], params[:np_z])
        else
            dim_error(params[:D])
        end
    end

    params = merge!(params, Dict(:procs_partition => procs_partition))

    if params[:case] == "Channel" && !params[:periodic] && params[:body_force] != 0.0
        println("For the non periodic channel is suggested not to have a body force")
    elseif params[:case] == "Channel" && params[:periodic] && params[:body_force] == 0.0
        println("For the  periodic channel is suggested to have a non-zero body force")
    end

    if params[:case] != "Channel" && params[:case] != "TaylorGreen"
        ν = params[:u_in] * params[:c] / params[:Re] #m2/s 
        params = merge!(params, Dict(:ν => ν))
        println("ν overwritten, computed from Reynolds $(params[:Re]), velocity $(params[:u_in]) and domain dimension $(params[:c]) new value: $(params[:ν])")
    end

    if params[:case] == "Channel"
        Lx = 2 * pi
        Ly = 2
        Lz = 2 / 3 * pi
        params = merge!(params, Dict(:Lx => Lx, :Ly => Ly, :Lz => Lz))
    end


    #Add parameters
    if params[:ode_method] == :AlphaMethod && !(:ρ∞ in keys(params))
        params = merge!(params, Dict(:ρ∞ => 0.8))

    elseif params[:ode_method] == :ThetaMethod && !(:θ in keys(params))
        params = merge!(params, Dict(:θ => 1))
    end


    if !(:nls_iter in keys(params))
        params = merge!(params, Dict(:nls_iter => 30)) #add a default value of nls iterations
    end

 
    if params[:solver]== :petsc && typeof(params[:options]) <: String
  
        if findfirst("snes", params[:options]) === nothing
            petsc_snes = false
        else
            petsc_snes = true
        end

        if findfirst("fieldsplit", params[:options]) === nothing
            petsc_fieldsplit = false
        else
            petsc_fieldsplit = true
        end

        params = merge!(params, Dict(:petsc_fieldsplit => petsc_fieldsplit))
        params = merge!(params, Dict(:petsc_snes => petsc_snes))

    end
    
    if params[:restart]
        restart_df = DataFrame(CSV.File(params[:restart_file]))
        params = merge!(params, Dict(:restart_df => restart_df))
    end

    if params[:linear]
        ũ_coeff = [2.1875, -2.1875, 1.3125, -0.3125]
        merge!(params, Dict(:ũ_coeff => ũ_coeff))
    end

    "Set the turbulence parameters, in specific the cache for generating the Eddies"
    if params[:case] =="Airfoil" || params[:case] =="Plate" || params[:case] =="Channel"
        Vboxinfo = params[:Vbox] #turbulence_box()
        Re = set_reynolds_stress(params::Dict)
        sem_cache = SEMcache(Vboxinfo, params[:dt], params[:u_in], Vboxinfo.σ, Re, params[:t0], nothing)
        params = merge!(params, Dict(:sem_cache => sem_cache))
    end
    
    
    
    
    
    if !(:debug_mode in keys(params))
        params = merge!(params, Dict(:debug_mode => false))
    elseif params[:debug_mode]
        if params[:solver] ==:petsc
            params[:options] = params[:options] * " -ksp_monitor -log_view"
        end
        params[:printmodel] = true
    end
    

    return params

end