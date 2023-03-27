#For periodic case look for ghost nodes, that's the isempty function; The y direction is not periodic for the channel; look for genralization?

"""

"""
function find_idx(p::VectorValue{2, Float64}, params; atol=1e-6)
  vv1 = findall(x->isapprox(x, p[1];  atol=atol), params[:restart_df].Points_0) 
  vv2 = findall(x->isapprox(x, p[2];  atol=atol), params[:restart_df].Points_1)

  if isempty(vv1) && params[:case] == "Channel"
    q = -sign(p[1]) * params[:Lx] + p[1]
    vv1 = findall(x-> isapprox(x, q;  atol=atol), params[:restart_df].Points_0) 
  end

  idx = vv1[findall(in(vv2),vv1)]
  if isempty(idx)
    idx = find_idx(p, params; atol = atol*10)
  end
  
  return idx[1]
end

function find_idx(p::VectorValue{3, Float64}, params; atol = 1e-4)
  vv1 = findall(x->isapprox(x, p[1];  atol=atol), params[:restart_df].Points_0) 
  vv2 = findall(x->isapprox(x, p[2];  atol=atol), params[:restart_df].Points_1)
  vv3 = findall(x->isapprox(x, p[3];  atol=atol), params[:restart_df].Points_2)
  
  if isempty(vv1) && params[:case] == "Channel"
    q = -sign(p[1]) * params[:Lx] + p[1]
    vv1 = findall(x-> isapprox(x, q;  atol=atol), params[:restart_df].Points_0) 
  end

  if isempty(vv3) && params[:case] == "Channel"
    q = -sign(p[3]) * params[:Lz] + p[3]
    vv3 = findall(x->isapprox(x, q; atol=atol) , params[:restart_df].Points_2) 
  end



  vv12 = vv1[findall(in(vv2),vv1)]
  idx = vv12[findall(in(vv3),vv12)]  
  
  if isempty(idx)
    idx = find_idx(p, params; atol = atol*10)
  end
  
  return idx[1]
end

function uh(p::VectorValue{2, Float64}, params::Dict{Symbol, Any}, idx::Int)

  VectorValue(params[:restart_df].uh_0[idx][1], params[:restart_df].uh_1[idx][1])
end


function uh(p::VectorValue{3, Float64}, params::Dict{Symbol, Any}, idx::Int)
  VectorValue(params[:restart_df].uh_0[idx][1], params[:restart_df].uh_1[idx][1], params[:restart_df].uh_2[idx][1])
end


function uh_restart(p, params::Dict{Symbol, Any})
  
  idx = find_idx(p, params)
  return uh(p, params, idx)

end

function ph_restart(p, params::Dict{Symbol, Any})
  idx = find_idx(p, params)
  ph = params[:restart_df].ph[idx][1]
  return ph
end


function restart_uh_field(params::Dict{Symbol, Any})
  println("u0")
  u0(x, t::Real) = uh_restart(x, params)
  u0(t::Real) = x -> u0(x, t::Real)
  return u0

end

function restart_ph_field(params::Dict{Symbol, Any})
  println("p0")
  
  p0(x, t::Real) =  ph_restart(x, params)
  p0(t::Real) = x -> p0(x, t::Real)

  return p0

end
