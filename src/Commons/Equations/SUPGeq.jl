
function h_param(Ω::GridapDistributed.DistributedTriangulation, D::Int64)
    h = map_parts(Ω.trians) do trian
        h_param(trian, D)
    end
    h = CellData.CellField(h, Ω)
    h
end


function h_param(Ω::Triangulation, D::Int64)
    h = lazy_map(h -> h^(1 / D), get_cell_measure(Ω))

    h
end

#Stabilization parameters
#Momentum stabilization
function τ(u, h, ν::Float64, dt::Float64)
        r = 1
        τ₂ = h^2 / (4 * ν)
        val(x) = x
        val(x::Gridap.Fields.ForwardDiff.Dual) = x.value
        u = val(norm(u))
    
        if iszero(u)
            return τ₂
        end
    
        τ₃ = dt / 2
        τ₁ = h / (2 * u)
        return 1 / (1 / τ₁^r + 1 / τ₂^r + 1 / τ₃^r)    
end
    
#Continuity stabilization
τb(u, h, ν::Float64, dt::Float64) = (u ⋅ u) * τ(u, h, ν, dt)



"""

SUPG non linear forumulation
"""
function SUPG(h, params::Dict{Symbol, Any})
    @unpack ν, dt, dΩ, hf, linear, θ = params
    if typeof(θ) <: AbstractVector
        θv,θp = θ
        θvp = θv/θp
    else
        θvp = 1
    end

    #Conservations equations
    #Conservations equations
    Rm, Rc = conservations_equations(params)
    dRm, dRc = derivative_conservations_equations(params)



    #Variational equations
    Bᴳ(t, (u, p), (v, q)) = ∫(∂t(u) ⋅ v)dΩ + ∫((u ⋅ ∇(u)) ⋅ v)dΩ - ∫((∇ ⋅ v) * p)dΩ + ∫((q * (∇ ⋅ u)))dΩ + ν * ∫(∇(v) ⊙ ∇(u))dΩ - ∫(hf(t) ⋅ v)dΩ

    #Stabilization equations
    B_stab(t, (u, p), (v, q)) = ∫((τ ∘ (u.cellfield, h, ν, dt) * (u ⋅ ∇(v) + ∇(q))) ⊙ Rm(t, (u, p)) # First term: SUPG, second term: PSPG u⋅∇(v) + ∇(q)
                                  +
                                  τb ∘ (u.cellfield, h, ν, dt) * (∇ ⋅ v) ⊙ Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
    )dΩ

    res(t, (u, p), (v, q)) = Bᴳ(t, (u, p), (v, q)) + B_stab(t, (u, p), (v, q))

    #SUPG Jacobian
    dBᴳ(t, (u, p), (du, dp), (v, q)) = ∫(((du ⋅ ∇(u)) ⋅ v) + ((u ⋅ ∇(du)) ⋅ v) + (∇(dp) ⋅ v) +  (q * (∇ ⋅ du)))dΩ + ν * ∫(∇(v) ⊙ ∇(du))dΩ
    dB_stab(t, (u, p), (du, dp), (v, q)) = ∫(((τ ∘ (u.cellfield, h, ν, dt) * (u ⋅ ∇(v)' +  ∇(q))) ⊙ dRm((u, p), (du, dp), (v, q))) + ((τ ∘ (u.cellfield, h, ν, dt) * (du ⋅ ∇(v)')) ⊙ Rm(t, (u, p))) + (τb ∘ (u.cellfield, h, ν, dt) * (∇ ⋅ v) ⊙ dRc(du)))dΩ


    jac(t, (u, p), (du, dp), (v, q)) = dBᴳ(t, (u, p), (du, dp), (v, q)) + dB_stab(t, (u, p), (du, dp), (v, q))

    #SUPG time-Jacobian                                   
    jac_t(t, (u, p), (dut, dpt), (v, q)) = ∫(dut ⋅ v)dΩ + ∫(τ ∘ (u.cellfield, h, ν, dt) * (u ⋅ ∇(v) + (θvp)*∇(q)) ⊙ dut)dΩ

    return res, jac, jac_t

    

end

"""

SUPG linearized
"""
function SUPG_lin(h, params::Dict{Symbol, Any})
    @unpack ν, dt, dΩ, hf, linear, ũ = params

    Rm, Rc = conservations_equations(params)

    a(t, (u, p), (v, q)) = ∫((ũ ⋅ ∇(u)) ⋅ v)dΩ - ∫((∇ ⋅ v) * p)dΩ + ∫((q * (∇ ⋅ u)))dΩ + ν * ∫(∇(v) ⊙ ∇(u))dΩ + 
    ∫((τ∘(ũ, h, ν, dt)*(ũ⋅(∇(v))' + ∇(q)))⊙Rm(u,p))dΩ +
    ∫(τb∘(ũ, h, ν, dt)*(∇⋅v)⊙Rc(u))dΩ 

    m(t, (u, p), (v, q)) =  ∫(u ⋅ v)dΩ + ∫((τ∘(ũ, h, ν, dt)*(ũ⋅(∇(v))' + ∇(q)))⊙(u))dΩ

    b(t, (v, q)) =∫(hf(t) ⋅ v)dΩ

    return m,a,b
end