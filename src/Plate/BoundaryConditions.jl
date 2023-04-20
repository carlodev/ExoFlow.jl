function bc_plate(params)

    @unpack u_in, D, xin, sem_cache = params
    
    
    # u_top(x, t) = (D == 2) ? VectorValue(u_in, 0.0) : VectorValue(u_in, 0, 0.0)
    # u_top(t::Real) = x -> u_top(x, t)
    u_top(x, t) = physical_blasius(x,xin+ x[1],params)
    u_top(t::Real) = x -> u_top(x, t)

    u_inlet(x, t) = physical_blasius(x,xin+ x[1],params) + generation_u_fluct!(x, t, sem_cache)
    u_inlet(t::Real) = x -> u_inlet(x, t)

    u_wall(x, t) = (D == 2) ? VectorValue(0.0, 0.0) : VectorValue(0.0, 0, 0.0)
    u_wall(t::Real) = x -> u_wall(x, t)

    u_start(x,t) = physical_blasius(x,xin + x[1],params)
    u_start(t::Real) = x -> u_start(x, t)

    if D == 2
        #Points
        #tag_1: bottom left
        #tag_2: bottom right
        #tag_3: top right
        #tag_4: top left
        #Lines
        #tag_5: horizonal bottom
        #tag_6: horizonal top
        #tag_7: veritcal left
        #tag_8: veritcal right
        u_diri_tags=["tag_1","tag_3", "tag_7", "tag_6", "tag_5", "tag_2", "tag_4"]
        u_diri_values = [u_top, u_top, u_inlet, u_top, u_wall, u_wall, u_top]
        p_diri_tags= ["tag_8"]
        p_diri_values = [0.0]
    elseif D == 3

        top = ["tag_20","tag_24"] # face

        bottom = ["tag_18", "tag_23"] # face
    
        inlet = ["tag_17", "tag_19", "tag_25"] #bottom, top, face
        outlet = ["tag_26"] #face
        

      
        u_diri_tags = append!(top, bottom, inlet)

        u_diri_values = [u_top, u_top, 
                        u_wall, u_wall,
                        u_top, u_top, u_top]

        p_diri_tags = outlet
      
        p_diri_values = [0.0]
        @assert length(u_diri_values) == length(u_diri_tags)
        @assert length(p_diri_values) == length(p_diri_tags)
        
        end

    return u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u_start

end