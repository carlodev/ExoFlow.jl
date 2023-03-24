#Functions for computing forces


function forces_domain(model, force_tags, degree)
    Γ = BoundaryTriangulation(model; tags=force_tags) 
    dΓ = Measure(Γ,degree)
    n_Γ = get_normal_vector(Γ) #Check if the direction is inside or outdside the domain; probably a -1 coefficient needed
    Γ, dΓ, n_Γ
end


function compute_forces(p, u, n_Γ, dΓ, ρ, ν, tn)
    μ = ρ*ν
    Fp = sum(∫(p⋅n_Γ)dΓ) #force from the pressure distribution
    Friction = μ*sum(∫(transpose(∇(u))⋅n_Γ)dΓ) #should be equivalent to ∂(u)/∂(n)
    #In a 3D case the forces in the z-direction are not computed 

    return [tn, Fp[1], Fp[2], Friction[1], Friction[2]]
end

function initialize_force_params(force_params, case, D, ν)
    if force_params !== nothing
    force_tags = force_params[:force_tags]
    ρ = force_params[:ρ]
    degree = force_params[:degree]
    model = force_params[:model]
    Γ, dΓ, n_Γ = forces_domain(model, force_tags, degree)  
    force_params = Dict(:Γ => Γ, :dΓ => dΓ, :n_Γ => n_Γ, :ρ => ρ,  :ν => ν, :output_file => "$(case)_$(D)d_forces.csv")      
    df = DataFrame(tn = 0.0, Fpx = 0.0, Fpy = 0.0, Frictionx = 0.0, Frictiony = 0.0 )
    CSV.write(force_params[:output_file], df, delim = ',')    
    
    #Initzialization
    return force_params

    else
        return nothing
    end

end


#Write force params
function write_forces(params::Dict{Symbol, Any}, uh_tn, ph_tn, tn::Float64)
    force_params = params[:force_params] #more compact

    if (force_params !== nothing) 
    Fout = compute_forces(ph_tn, uh_tn, force_params[:n_Γ], force_params[:dΓ], params[:ρ], params[:ν], tn)
    println(Fout)
    df = DataFrame(tn=Fout[1], Fpx=Fout[2], Fpy=Fout[3], Frictionx=Fout[4], Frictiony=Fout[5])
    CSV.write(force_params[:output_file], df, delim=',', append=true)
end

end