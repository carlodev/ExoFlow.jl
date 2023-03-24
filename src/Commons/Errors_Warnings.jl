function dim_error(D)
    error("Dimension $D not supported, only 2 or 3 dimensional cases")
end



function procs_error(D, nn)

    if length(nn) != D
        error("Partition non valid")

    end

end

function method_error(method)
    if (method != :VMS) && (method != :SUPG)
        error("Method $method not supported")
    end
end

function solver_error(solver, options)
    if (solver == :petsc) && (options == "")
        error("Specify options for petsc solver")
    end
end


function dim_not_supported(D, case)

    if D != 2
        error("$case dimension D = $D not supported, only 2-D is supported")
    end
end



#Warnings

function params_integrity(params)
    method_error(params[:method])
    solver_error(params[:solver], params[:options])
end

