include("DNS_start.jl")

function initial_velocity(params)
    @unpack D, start_condition, sem_cache = params
    @unpack u_in, periodic, ν, body_force = params

    #Poiseuille parabolic profile, for h = 1 
    umax = 1 ./ (2 .* ν) .* body_force
    ux_Poiseuille(x) = umax .* (1 .- x[2] .^ 2)
    ux_laminar(x) = (D == 2) ? VectorValue(ux_Poiseuille(x), 0.0) : VectorValue(ux_Poiseuille(x), 0.0, 0.0) #Because h=1, simplified expression

    #Turbulent Average Profile from DNS. It gives a VectorValue results
    u_turbulent = dns_velocity()

    #Chosing the average profile
    u_profile(x) = (start_condition == :laminar) ? ux_laminar(x)  : u_turbulent(x)


    ux_noise(x) = (rand(1)[1] .- 0.5) .* 0.3 .* u_profile(x)[1] #Or just .* 2/3.* umax #more turbulent/chaotic
    uy_noise(x) = (rand(1)[1] .- 0.5) .* 0.2 #are function of x, in this way at each point there is a different value, otherwise it is computed just once
    uz_noise(x) = (rand(1)[1] .- 0.5) .* 0.2

    u_noise(x, t::Real) = (D == 2) ? VectorValue(ux_noise(x), uy_noise(x)) : VectorValue(ux_noise(x), uy_noise(x), uz_noise(x))
    u0_turb(x, t::Real) = u_profile(x) + u_noise(x,t) #generation_u_fluct!(x, t, sem_cache)

    u0(x, t::Real) = periodic ? u0_turb(x, t) : ((D == 2) ? VectorValue(u_in, 0.0) : VectorValue(u_in, 0.0, 0.0))
    u0(t::Real) = x -> u0(x, t::Real)

    return u0

end





function bc_channel(params)


    u_walls(x, t::Real) = (params[:D] == 2) ? VectorValue(0.0, 0.0) : VectorValue(0.0, 0.0, 0.0)
    u_walls(t::Real) = x -> u_walls(x, t)

    p_diri_tags = String[]
    p_diri_values = Float64[]

    if params[:periodic] #in periodic case pressure zero fixed in the centre point of the domain
        p_diri_tags = "centre"
        p_diri_values = 0.0
    end

    if params[:D] == 2
        top = "tag_5"
        bottom = "tag_6"
        inlet = "tag_7"
        outlet = "tag_8"
        outlet_top = "tag_2" # top right corner, not adding corners results in a "jump" at the outlet
        outlet_bottom = "tag_4" # bottom right corner
        inlet_top = "tag_1"
        inlet_bottom = "tag_3"
        u_diri_tags = [top, bottom]



        if params[:periodic]
            u_diri_values = [u_walls, u_walls]

        else

            append!(u_diri_tags, [inlet, inlet_top, inlet_bottom, outlet_top, outlet_bottom])
            append!(p_diri_tags, [outlet])
            u_diri_values = [u_walls, u_walls, params[:u0], u_walls, u_walls, u_walls, u_walls]
            append!(p_diri_values, [0.0])
        end



    elseif params[:D] == 3
        top = ["tag_23"] # face
        top_sides = ["tag_09", "tag_11"] # right left -not in periodic boundaries
        bottom = ["tag_24"] # face
        bottom_sides = ["tag_10", "tag_12"]

        inlet = ["tag_25"]
        inlet_corners = ["tag_01", "tag_03", "tag_05", "tag_07"] #topleft bottomleft topright bottomright
        inlet_sides = ["tag_13", "tag_15", "tag_17", "tag_19"] #left right top bottom

        outlet = ["tag_26"]
        outlet_corners = ["tag_02", "tag_04", "tag_06", "tag_08"] #topleft bottomleft topright bottomright
        outlet_sides = ["tag_14", "tag_16", "tag_18", "tag_20"]

        sides = ["tag_21", "tag_22"] # left right
        u_diri_tags = append!(top, bottom)



        if params[:periodic]
            u_diri_values = [u_walls, u_walls]

        else
            append!(u_diri_tags, top_sides, bottom_sides, inlet, inlet_corners, inlet_sides, outlet_corners, outlet_sides, sides)
            append!(p_diri_tags, outlet)
            u_diri_values = [u_walls, u_walls, u_walls, u_walls, u_walls, u_walls,
                params[:u0],
                params[:u0], params[:u0], params[:u0], params[:u0],
                params[:u0], params[:u0], params[:u0], params[:u0],
                u_walls, u_walls, u_walls, u_walls,
                u_walls, u_walls, u_walls, u_walls,
                u_walls, u_walls]

            append!(p_diri_values, 0)
        end
    else
        dim_error(params[:D])
    end

    return u_diri_tags, u_diri_values, p_diri_tags, p_diri_values


end
