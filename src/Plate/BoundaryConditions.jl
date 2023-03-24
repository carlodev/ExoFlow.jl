function bc_plate(params)

    u_in = params[:u_in]
    
    
    u_top(x, t) = VectorValue(u_in, 0)
    u_top(t::Real) = x -> u_top(x, t)

    u_wall(x, t) = (x[1]<0.5) ? u_top(x, t) : VectorValue(0,0)
    u_wall(t::Real) = x -> u_wall(x, t)



    
    u_diri_tags=["tag_1","tag_3", "tag_7", "tag_6", "tag_5", "tag_2"]
    u_diri_values = [u_top, u_top, u_top, u_top, u_wall, u_wall]
    p_diri_tags= ["tag_8", "tag_4"]
    p_diri_values = [0.0, 0.0]
   

    return u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u_top

end