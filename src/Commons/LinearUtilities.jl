#MPI Backend
"""
It instantiate an AbstactVector of 4 elements. Each element is an AbstractVector where the values of the velocity in each node is stored.
The first element refers to the time step n, the second to the time step n-1, the third n-2 and the last n-3.
"""
function create_ũ_vector(zfields::MPIData)
    zfv = deepcopy(get_part(zfields))
    zfv1 = get_free_dof_values(zfv)
    return [zfv1, zfv1, zfv1, zfv1]
end

"""
It updates the vector which stores the values of velocity at previous time steps.
"""
function  update_ũ_vector!(ũ_vec::Vector{Vector{Float64}}, uhfields::MPIData)
    uh_new = get_free_dof_values(get_part(uhfields))
    circshift!(ũ_vec,1)
    ũ_vec[1] = deepcopy(uh_new)
end

"""
It updates the convective velocity exitmation ``\\tilde{u}`` for the approximation of the non linear term: ``\\tilde{u} \\cdot \\nabla(u)``
"""
function update_ũ(ũ_vec::Vector{Vector{Float64}}, coeff::Vector{Float64})
    updt_ũ = ũ_vec[1]*coeff[1] + ũ_vec[2] *coeff[2] + ũ_vec[3] *coeff[3] + ũ_vec[4]*coeff[4]
    return updt_ũ
end


function update_free_values!(zfields::MPIData, zt::Vector{Float64})
    copyto!(zfields.part.free_values, zt)
end




#Sequential Backend
"""
  create_ũ_vector(zfields::SequentialData)
"""
function create_ũ_vector(zfields::SequentialData)
    u_vec = Vector[]
    for p = 1:1:length(zfields.parts)
      zfv = get_free_dof_values(get_part(zfields,p))
      push!(u_vec, [zfv, zfv, zfv, zfv])
    end
  return u_vec  
end

function update_ũ_vector!(ũ_vec::Vector{Vector}, zfields::SequentialData)
for p = 1:1:length(zfields.parts)
      zfv = get_free_dof_values(get_part(zfields,p))
      circshift!(ũ_vec[p],1)
      ũ_vec[p][1] = deepcopy(zfv)
    end
end
  

function update_ũ(ũ_vec::Vector{Vector}, coeff::Vector{Float64})
    updt_ũ = Vector[]
    for u_vec in ũ_vec
    updt_u = u_vec[1]*coeff[1] + u_vec[2] *coeff[2] + u_vec[3] *coeff[3] + u_vec[4]*coeff[4]
    push!(updt_ũ, updt_u)
    end
  
    return updt_ũ
end


function update_free_values!(zfields::SequentialData, zt::Vector{Vector})
    for p = 1:1:length(zfields.parts)
      copyto!(zfields.parts[p].free_values, zt[p])
    end
end


"""
    update_linear!(params::Dict{Symbol,Any})

Wrapper function for updating the ũ_vector containing the speed values at previous time steps 
and ũ which is the approximation of ũ for the non linear terms, ũ⋅(∇u)
"""
function update_linear!(params::Dict{Symbol,Any}, uh_tn)
  if params[:linear]
    update_ũ_vector!(params[:ũ_vector], uh_tn.fields)
    zt = update_ũ(params[:ũ_vector], params[:ũ_coeff])
    update_free_values!(params[:ũ].fields, zt)
  end
end


#REPL
using Gridap.FESpaces
"""
  create_ũ_vector(zfields::SingleFieldFEFunction)
"""
function create_ũ_vector(zfields::SingleFieldFEFunction)
    u_vec = Vector[]
    zfv = get_free_dof_values(zfields)
    u_vec =  [zfv, zfv, zfv, zfv]
  return u_vec  
end

function update_ũ_vector!(ũ_vec::Vector{Vector{Float64}}, zfields::SingleFieldFEFunction)
      zfv = get_free_dof_values(zfields)
      circshift!(ũ_vec,1)
      ũ_vec[1] = deepcopy(zfv)
end
  


function update_free_values!(zfields::SingleFieldFEFunction, zt::Vector{Float64})
  copyto!(zfields.free_values, zt)
end

function update_linear!(params::Dict{Symbol,Any}, uh_tn::SingleFieldFEFunction)
  if params[:linear]
    update_ũ_vector!(params[:ũ_vector], uh_tn)
    zt = update_ũ(params[:ũ_vector], params[:ũ_coeff])
    update_free_values!(params[:ũ], zt)
  end
end

