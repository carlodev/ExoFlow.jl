function bc_liddriven(model, params)
    u_in = params[:u_in]
    D = params[:D]
    mesh_gen = params[:mesh_gen]
    labels = get_face_labeling(model)
    
    u_wall(x,t) = (D == 2) ? VectorValue(0.0, 0.0) : VectorValue(0.0, 0.0, 0.0) 
    u_wall(t::Real) = x -> u_wall(x,t)

    
    u_top(x,t) = (D == 2) ? VectorValue(u_in, 0.0) : VectorValue(u_in, 0.0, 0.0) 
    u_top(t::Real) = x -> u_top(x,t)

    if mesh_gen
        
        u_diri_tags=["top", "walls", "walls_p"]
        u_diri_values = [u_top, u_wall, u_wall]
        p_diri_tags=  "walls_p"

    elseif D == 2
        add_tag_from_tags!(labels, "diri1", [5,])
        add_tag_from_tags!(labels, "diri0", [1, 2, 4, 3, 6, 7, 8])
        add_tag_from_tags!(labels, "p", [4,])
        u_diri_tags=["diri0", "diri1"]
        u_diri_values = [u_wall, u_top]
        p_diri_tags= "p"
    
    elseif D == 3
        println("Model 3 dimension")

        add_tag_from_tags!(labels, "diri1", [3,4,7,8,10,12,19,20,24])
        add_tag_from_tags!(labels, "diri0", [1, 2,5,6,9,11,13,14,15,16,17,18,21,22,23,25,26])
        add_tag_from_tags!(labels, "p", [1,])
        u_diri_tags=["diri0", "diri1"]
        u_diri_values = [u_wall, u_top]
        p_diri_tags= "p"

    end
  

    p_diri_values = 0.0

    return u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u_wall

end