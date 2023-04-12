"""
    conservations_equations(params::Dict{Symbol,Any})

It provides the Navier Stokes equations: continuity and momentum. The equations are different if the problem is linearized or not.
Conitnuity:
``\\nabla\\cdot u = 0``

Momentum:
``\\dfrac{\\partial u}{\\partial t} + u \\cdot \\nabla(u) + \\nabla (p) + \\nu \\Delta u= 0``

Momentum linearized:
``\\tilde{u} \\cdot \\nabla(u) + \\nabla (p) + \\nu \\Delta u = 0``
where ``\\tilde{u}`` is the approximation of the convective velocity.
"""
function conservations_equations(params::Dict{Symbol,Any})
    #Momentum
    @unpack linear, steady, ν, hf= params
    if linear
        @unpack ũ = params
    end

    println("linear = $linear, steady = $steady")

    #Linearized momentum equation
    # Momentum residual, without viscous term if first order elements
    skew_lin(u) = ũ *(∇⋅u)/2
    Rm_lin(u, p) = ũ ⋅ ∇(u) + ∇(p) + skew_lin(u)#- ν*Δ(u)       
    
    # Momentum residual, without viscous term if first order elements
    # Rm_steady(u, p) = u ⋅ ∇(u) + ∇(p) - hf - ν*Δ(u)

    Rm(t, (u, p)) = ∂t(u) + u ⋅ ∇(u) + ∇(p) - hf(t)  #- ν*Δ(u)

 


    # Continuity residual
    Rc(u) = ∇ ⋅ u

    if linear
        return Rm_lin, Rc
    #elseif steady
    #   return Rm_steady, Rc
    else
        return  Rm,  Rc
    end

end

"""
    derivative_conservations_equations(params::Dict{Symbol,Any})

It provides the spatial derivative of the Navier Stokes equations: continuity and momentum. The equations are different if the problem is linearized or not.
Conitnuity:
``\\nabla\\cdot du = 0``

Momentum:
`` du \\cdot \\nabla(u) + u \\cdot \\nabla(du) + \\nabla (dp) + \\nu \\Delta du= 0``

Momentum linearized:
``\\tilde{u} \\cdot \\nabla(u) + \\nabla (p) + \\nu \\Delta u = 0``
"""
function derivative_conservations_equations(params::Dict{Symbol,Any})
    @unpack linear, ν = params
    if linear
        @unpack ũ = params
    end
    #Linearized momentum equation
    #Momentum derivative
    dRm_lin((u, p), (du, dp), (v, q)) = ũ ⋅ ∇(du) + ∇(dp) #- ν*Δ(du)

    #Momentum derivative
    dRm((u, p), (du, dp), (v, q)) = du ⋅ ∇(u) + u ⋅ ∇(du) + ∇(dp) #- ν*Δ(du)

    # Continuity derivative residual

    dRc(du) = ∇ ⋅ du

    if linear
        return dRm_lin, dRc
    else
        return dRm, dRc
    end

end


