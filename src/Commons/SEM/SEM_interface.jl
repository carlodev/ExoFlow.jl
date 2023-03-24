mutable struct SEMcache
    Vboxinfo::VirtualBox
    dt::Float64
    U₀::Float64
    σ::Vector{Float64}
    Re::Union{Matrix,Reynolds_stress_interpolator}
    t_old::Float64
    Eddies::Union{Vector{SemEddy},Nothing}
end

function generation_u_fluct!(x::VectorValue, t::Real, sem_cache::SEMcache)
    
    if sem_cache.Eddies === nothing 
        sem_cache.Eddies = initialize_eddies(sem_cache.Vboxinfo)
    end
    
    
    if t>sem_cache.t_old
        #convecting eddies
        @info "convecting eddies"
        sem_cache.Eddies = map(ej -> convect_eddy(sem_cache.dt, ej, sem_cache.U₀, sem_cache.Vboxinfo), sem_cache.Eddies)
        # comm = MPI.COMM_WORLD
        # MPI.Barrier(comm)
        #Update cache: time and Eddies
        sem_cache.t_old = t
    
    elseif t<sem_cache.t_old
        error("cannot compute eddies at a previous time step")
    end
    

    D = length(x)

    if D ==3
        vec_point = [x[1], abs(x[2]), x[3]]

    elseif D == 2
        vec_point = [x[1], abs(x[2]), 0.0]

    end
    Re_point = SyntheticEddyMethod.Reynolds_stress_points(vec_point, sem_cache.Re)
    
    #compute fluctuation. U₀ = 0.0, added in the boundary condtions
    @time U = compute_uDFSEM(vec_point, sem_cache.Eddies, sem_cache.Vboxinfo, Re_point)[1]
        
    if D == 3
        return VectorValue(U...)
    
    elseif D== 2
        return VectorValue(U[1], U[2])
    
    end

end

function create_Re_stress(params::Dict{Symbol, Any})
@unpack Re_filename, TI = params
if Re_filename == "none"
    @unpack u_in = params
    u_p = (u_in * TI)^2
    Re = collect(I(3).*u_p)
 else
    @unpack dims = params
    Re = get_reynolds_stress_from_file(Re_filename; dims = dims)
end

return Re
end

