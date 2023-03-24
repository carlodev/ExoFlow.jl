function τmh(uu, G, GG)
    println("inside tmh0")

    τ₁ =  (2 *Cᵢ[1] / dt)^2 #Here, you can increse the 2 if CFL high
    τ₃ = Cᵢ[2] * (ν^2 * GG)
    println("inside tmh")
    val(x) = x
    function val(x::Gridap.Fields.ForwardDiff.Dual)
      x.value
    end
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
    println("inside tmh 1")
  
    if iszero(norm(uu_new))
      return (τ₁ .+ τ₃) .^ (-1 / 2)
    end
    println("inside tmh 2")
  
    τ₂ = uu_new ⋅ G ⋅ uu_new
    return (τ₁ .+ τ₂ .+ τ₃) .^ (-1 / 2)
  end
  
  function τch(uu, gg, G, GG)
    return 1 / (τmh(uu, G, GG) ⋅ gg)
  end

  """
  if params[:method] == :VMS
      @unpack ν, dt, Cᵢ, D, G, GG, gg = params
      function τmh(uu, G, GG)
          println("inside tmh0")
      
          τ₁ =  (2 *Cᵢ[1] / dt)^2 #Here, you can increse the 2 if CFL high
          τ₃ = Cᵢ[2] * (ν^2 * GG)
          println("inside tmh")
          val(x) = x
          function val(x::Gridap.Fields.ForwardDiff.Dual)
            x.value
          end
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
          println("inside tmh 1")
        
          if iszero(norm(uu_new))
            return (τ₁ .+ τ₃) .^ (-1 / 2)
          end
          println("inside tmh 2")
        
          τ₂ = uu_new ⋅ G ⋅ uu_new
          return (τ₁ .+ τ₂ .+ τ₃) .^ (-1 / 2)
        end
        
        function τch(uu, gg, G, GG)
          return 1 / (τmh(uu, G, GG) ⋅ gg)
        end
      
        



      println("Computing stabilization parameters")
      τm_h = τmh∘(uh_tn, G, GG)
      println("tm")
      τc_h = τch∘(uh_tn,  gg, G, GG)
      println("tc")

      pvd[tn] = createvtk(params[:Ω], joinpath("Results", "$(params[:case])_$(params[:D])d_$tn" * ".vtu"), cellfields=["uh" => uh_tn, "ph" => ph_tn, "wh" => ωh_tn, "tau_m"=>τm_h, "tau_c"=>τc_h])
  end
  """