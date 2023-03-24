function conservations_equations(params)
    #Momentum
    @unpack linear, steady, ν, hf= params
    if linear
        @unpack ũ = params
    end

    println("linear = $linear, steady = $steady")

    #Linearized momentum equation
    # Momentum residual, without viscous term if first order elements
    Rm_lin(u, p) = ũ ⋅ ∇(u) + ∇(p) #- ν*Δ(u)        

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


function derivative_conservations_equations(params)
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




function time_derivative_conservations_equation()
    

end


