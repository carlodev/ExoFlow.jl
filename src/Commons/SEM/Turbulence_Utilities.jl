
function set_reynolds_stress(params::Dict{Symbol,Any})

@unpack Re_filename, TI = params
if Re_filename == "none"
    @unpack u_in = params
    Re = set_reynolds_stress(u_in, TI)
 else
    @unpack dims = params
    Re = set_reynolds_stress(Re_filename, dims)
end

return Re

end



function set_reynolds_stress(U₀::Float64, TI::Float64)
    #set diagonal reynolds stress the parameters
    u_p = (U₀ * TI)^2
    Re_stress = collect(I(3).*u_p)

    return  Re_stress
    
end

function set_reynolds_stress(reynolds_stress_file::String, dims::Tuple)
    #set diagonal reynolds stress from file
    Re_stress = get_reynolds_stress_from_file(reynolds_stress_file; dims = dims)

    return  Re_stress
    
end


function set_reynolds_stress()
    #Manually set the reynolds stress
    Re_stress = rand(3,3)

    return  Re_stress
end