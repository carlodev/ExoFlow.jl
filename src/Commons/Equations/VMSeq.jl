
"""
VMS implementation, based on Bazilevs et al. DOI: 10.1016/j.cma.2007.07.016
"""

#VMS parameters

"""
It provides the set of stabilization parameters for the VMS: G, GG, gg. 
GG comes from the inverse of the jacobian of the map between reference and physical domain. It is exactly the cell dimension obtained by [`h_param`](@ref) for orthogonal carthesian grids.

"""
function G_params(Ω::GridapDistributed.DistributedTriangulation, params)
  @time G, GG, gg = map_parts(Ω.trians) do trian

    G_params(trian, params)
  end

  G = CellData.CellField(G, Ω)
  GG = CellData.CellField(GG, Ω)
  gg = CellData.CellField(gg, Ω)
  G, GG, gg
end

function G_params(trian::Gridap.Geometry.BodyFittedTriangulation, params) #trian == Ω
  D = params[:D]

  ξₖ = get_cell_map(trian)
  Jt = lazy_map(Broadcasting(∇), ξₖ)
  inv_Jt = lazy_map(Operation(inv), Jt)

  if D == 2
    eval_point = Point(0.5, 0.5)
  else
    eval_point = Point(0.5, 0.5, 0.5)
  end

  d = lazy_map(evaluate, inv_Jt, Fill(eval_point, num_cells(trian)))
  dtrans = lazy_map(Broadcasting(transpose), d)
  G = lazy_map(Broadcasting(⋅), d, dtrans)
  GG = lazy_map(Broadcasting(⊙), G, G)


  function gg_operation(d)

    if D == 2

      return (d[1] + d[3])^2 + (d[2] + d[4])^2


    elseif D == 3

      return (d[1] + d[4] + d[7])^2 + (d[2] + d[5] + d[8])^2 + (d[3] + d[6] + d[9])^2

    end

  end

  gg = lazy_map(Broadcasting(gg_operation), d)

  print2file_stab_params(params, [G, GG, gg])

  G, GG, gg
end #end G_params - Triangulation


#Stabilization paramameter momentum
"""
Stabilization parameter for momentum equation for the VMS.
"""
  function τm(uu, G, GG, Cᵢ1::Int64, Cᵢ2::Int64,  ν::Float64,dt::Float64)
    τ₁ = Cᵢ1 * (2 / dt)^2 #Here, you can increse the 2 if CFL high
    τ₃ = Cᵢ2 * (ν^2 * GG)

    val(x) = x
    function val(x::Gridap.Fields.ForwardDiff.Dual)
      x.value
    end
    D = length(uu)
    if D == 2
      uu1 = val(uu[1])
      uu2 = val(uu[2])
      uu_new = VectorValue(uu1, uu2)
    elseif D == 3
      uu1 = val(uu[1])
      uu2 = val(uu[2])
      uu3 = val(uu[3])
      uu_new = VectorValue(uu1, uu2, uu3)
    end

    if iszero(norm(uu_new))
      return (τ₁ .+ τ₃) .^ (-1 / 2)
    end

    τ₂ = uu_new ⋅ G ⋅ uu_new
    return (τ₁ .+ τ₂ .+ τ₃) .^ (-1 / 2)
  end

#Stabilization paramameter continuity
"""
Stabilization parameter for contitnuity equation for the VMS.
"""
function τc(uu, gg, G, GG, Cᵢ1::Int64, Cᵢ2::Int64, ν,dt)
    return 1 / (τm(uu, G, GG,  Cᵢ1, Cᵢ2, ν,dt) ⋅ gg)
end

"""
VMS non linear variational forumulation. It calles the [`conservations_equations`](@ref) and [`derivative_conservations_equations`](@ref).
It provides the equations set, the jacobian and the time jacobian.
"""
function VMS(G, GG, gg, params::Dict{Symbol,Any})
  @unpack ν, dt, dΩ, hf, D, Cᵢ, θ = params

  if typeof(θ) <: AbstractVector
    θv,θp = θ
    θvp = θv/θp
else
    θvp = 1
end

Cᵢ1 = Cᵢ[1]
Cᵢ2 = Cᵢ[2]

  if params[:debug_mode]
    merge!(params, Dict(:G => G, :GG => GG, :gg => gg))
  end

  #Conservations equations
  Rm, Rc = conservations_equations(params)
  dRm, dRc = derivative_conservations_equations(params)

  dtRm((u, p), (dut, dpt), (v, q)) = dut 



  #VMS equations
  #Function used in VMS terms
  TRm(t, (u, p)) = τm ∘ (u.cellfield, G, GG,  Cᵢ1, Cᵢ2, ν,dt) * Rm(t, (u, p))

  #Variational equations
  Bᴳ(t, (u, p), (v, q)) = ∫(∂t(u) ⋅ v)dΩ + ∫((u ⋅ ∇(u)) ⋅ v)dΩ - ∫((∇ ⋅ v) * p)dΩ + ∫((q * (∇ ⋅ u)))dΩ + ν * ∫(∇(v) ⊙ ∇(u))dΩ - ∫(hf(t) ⋅ v)dΩ

  #SUPG terms 
  B_SUPG(t, (u, p), (v, q)) = ∫((u ⋅ ∇(v) + ∇(q)) ⊙ TRm(t, (u, p)))dΩ + ∫((∇ ⋅ v) ⊙ (τc ∘ (u.cellfield, gg, G, GG,  Cᵢ1, Cᵢ2, ν,dt) * Rc(u)))dΩ

  #First VMS term
  B_VMS1(t, (u, p), (v, q)) = ∫((u ⋅ (∇(v))') ⊙ TRm(t, (u, p)))dΩ

  #Second VMS term
  B_VMS2(t, (u, p), (v, q)) = -1 * ∫((∇(v)) ⊙ (outer(TRm(t, (u, p)), TRm(t, (u, p)))))dΩ

  #Adding all the contributions
  Bᴹ(t, (u, p), (v, q), ν, dΩ) = Bᴳ(t, (u, p), (v, q)) + B_SUPG(t, (u, p), (v, q)) + B_VMS1(t, (u, p), (v, q)) + B_VMS2(t, (u, p), (v, q))

  res(t, (u, p), (v, q)) = Bᴹ(t, (u, p), (v, q), ν, dΩ)

  #VMS Jacobian
  #Function used in the VMS terms -  derivative
  dTRm((u, p), (du, dp), (v, q)) = τm ∘ (u.cellfield, G, GG, Cᵢ1, Cᵢ2, ν,dt) * dRm((u, p), (du, dp), (v, q)) #Derivative of the Function used in the second VMS term

  #Variational equations - derivative
  dBᴳ(t, (u, p), (du, dp), (v, q)) = ∫(((du ⋅ ∇(u)) ⋅ v) + ((u ⋅ ∇(du)) ⋅ v) + (∇(dp) ⋅ v) + (q * (∇ ⋅ du)))dΩ + ν * ∫(∇(v) ⊙ ∇(du))dΩ

  #SUPG terms derivative
  dB_SUPG(t, (u, p), (du, dp), (v, q)) = ∫(((du ⋅ ∇(v)) ⋅ TRm(t, (u, p))) + ((u ⋅ ∇(v) + ∇(q)) ⋅ dTRm((u, p), (du, dp), (v, q))) + ((∇ ⋅ v) ⋅ (τc ∘ (u.cellfield, gg, G, GG, Cᵢ1, Cᵢ2, ν,dt) .* dRc(du))))dΩ

  #First VMS term derivative
  dB_VMS1(t, (u, p), (du, dp), (v, q)) = ∫(((du ⋅ ∇(v)') ⊙ TRm(t, (u, p))) + ((u ⋅ ∇(v)') ⊙ dTRm((u, p), (du, dp), (v, q))))dΩ


  #Second VMS term derivative
  dB_VMS2(t, (u, p), (du, dp), (v, q)) = ∫(∇(v) ⊙ (outer(dTRm((u, p), (du, dp), (v, q)), TRm(t, (u, p)))))dΩ + ∫(∇(v) ⊙ (outer(TRm(t, (u, p)), dTRm((u, p), (du, dp), (v, q)))))dΩ


  #Adding all the contributions
  jac(t, (u, p), (du, dp), (v, q)) = dBᴳ(t, (u, p), (du, dp), (v, q)) + dB_SUPG(t, (u, p), (du, dp), (v, q)) + dB_VMS1(t, (u, p), (du, dp), (v, q)) - dB_VMS2(t, (u, p), (du, dp), (v, q))

  #VMS time-Jacobian
  dtTRm(t, (u, p), (dut, dpt), (v, q)) = (τm ∘ (u.cellfield, G, GG, Cᵢ1, Cᵢ2, ν,dt)) ⋅ dut #Derivative of the Function used in the second VMS term

  dtBᴳ(t, (u, p), (dut, dpt), (v, q)) = ∫(dut ⋅ v)dΩ 
  dtB_SUPG(t, (u, p), (dut, dpt), (v, q)) =  ∫((u ⋅ ∇(v) + (θvp)*∇(q)) ⊙ dtTRm(t, (u, p), (dut, dpt), (v, q)))dΩ 
  dtB_VMS1(t, (u, p), (dut, dpt), (v, q)) = ∫((u ⋅ (∇(v))') ⊙ dtTRm(t, (u, p), (dut, dpt), (v, q)))dΩ
  dtB_VMS2(t, (u, p), (dut, dpt), (v, q)) =  ∫((∇(v)) ⊙ (outer(dtTRm(t, (u, p), (dut, dpt), (v, q)), TRm(t, (u, p)))))dΩ + ∫((∇(v)) ⊙ (outer(TRm(t, (u, p)), dtTRm(t, (u, p), (dut, dpt), (v, q)))))dΩ        
  
  jac_t(t, (u, p), (dut, dpt), (v, q)) = dtBᴳ(t, (u, p), (dut, dpt), (v, q)) +  dtB_SUPG(t, (u, p), (dut, dpt), (v, q)) + dtB_VMS1(t, (u, p), (dut, dpt), (v, q)) - dtB_VMS2(t, (u, p), (dut, dpt), (v, q))

  res, jac, jac_t

end #end VMS


"""
VMS linear variational forumulation. It calles the [`conservations_equations`](@ref) and [`derivative_conservations_equations`](@ref). 
The terms are divided between mass matrix, stifness matrix and right hand side.
"""
function VMS_lin(G, GG, gg, params::Dict{Symbol,Any})
  @unpack ν, dt, dΩ, hf, linear, ũ, Cᵢ, θ = params
  Cᵢ1 = Cᵢ[1]
  Cᵢ2 = Cᵢ[2]
  if typeof(θ) <: AbstractVector
    θv,θp = θ
    θvp = θv/θp
  else
    θvp = 1
  end

  aRm, aRc = conservations_equations(params)
  aTRm(u, p) = τm ∘ (ũ, G, GG,  Cᵢ1, Cᵢ2, ν,dt) * aRm(u, p)
  bTRm(u, p) = τm ∘ (ũ, G, GG,  Cᵢ1, Cᵢ2, ν, dt) * u
  bbTRm(u, p) = τm ∘ (ũ, G, GG,  Cᵢ1, Cᵢ2, ν,dt) * ũ
  
  fTRm(t) = τm ∘ (ũ, G, GG,  Cᵢ1, Cᵢ2, ν,dt) * hf(t)

  m(t, (u, p), (v, q)) =  ∫(u ⋅ v)dΩ + 
  ∫((ũ ⋅ ∇(v) + (θvp)*∇(q)) ⊙ bTRm(u, p))dΩ + 
  ∫((ũ ⋅ (∇(v))') ⊙ bTRm(u, p))dΩ #+
  #  -0.5 * ∫((∇(v)) ⊙ (outer(bTRm(u, p), bbTRm(u, p))))dΩ +
  #  -0.5 * ∫((∇(v)) ⊙ (outer(bbTRm(u, p), bTRm(u, p))))dΩ
  
  
  aᴳ(t, (u, p), (v, q)) = ∫((ũ ⋅ ∇(u)) ⋅ v)dΩ - ∫((∇ ⋅ v) * p)dΩ + ∫((q * (∇ ⋅ u)))dΩ + ν * ∫(∇(v) ⊙ ∇(u))dΩ 
  
  a_SUPG(t, (u, p), (v, q)) = ∫((ũ ⋅ ∇(v) + ∇(q)) ⊙ aTRm(u, p))dΩ + ∫((∇ ⋅ v) ⊙ (τc ∘ (ũ, gg, G, GG,  Cᵢ1, Cᵢ2, ν,dt) * aRc(u)))dΩ
  a_VMS1(t,(u, p), (v, q)) = ∫((ũ ⋅ (∇(v))') ⊙ aTRm(u, p))dΩ
  a_VMS2(t, (u, p), (v, q)) = -0.5 * ∫((∇(v)) ⊙ (outer(aTRm(u, p), bbTRm(u, p))))dΩ -0.5* ∫((∇(v)) ⊙ (outer(bbTRm(u, p), aTRm(u, p))))dΩ
  
  a(t, (u, p), (v, q)) = aᴳ(t, (u, p), (v, q)) +  a_SUPG(t, (u, p), (v, q)) + a_VMS1(t,(u, p), (v, q)) #+ a_VMS2(t,(u, p), (v, q))
  b(t, (v, q)) =∫(hf(t) ⋅ v)dΩ + ∫((ũ ⋅ ∇(v) + ∇(q)) ⋅  fTRm(t))dΩ + ∫((ũ ⋅ (∇(v))') ⋅ fTRm(t))dΩ

  return m,a,b
end

