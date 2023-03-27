"""
Add the body force. The user in the initial parameters specify the component in the x direction. Then the body force function is created 2D or 3D accordingly to the physics of the case.
"""
function hf_gen!(params::Dict{Symbol,Any})
    hf(x, t::Union{Real,AbstractVector}) = (params[:D] == 2) ? VectorValue(params[:body_force], 0.0) : VectorValue(params[:body_force], 0.0, 0.0)
    hf(t::Union{Real,AbstractVector}) = x -> hf(x, t)
    merge!(params, Dict(:hf => hf))
    return params
end

"""
It prints the mesh of the model
"""
function printmodel(params::Dict{Symbol,Any}, model)
    if params[:printmodel]
        @info "printing model"
        writevtk(model, "$(params[:case])_$(params[:D])d")
    end
end



"""
It creates the finite elements spaces accordingly to the previously generated dirichelet tags
"""
function creation_fe_spaces(params::Dict{Symbol,Any}, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{params[:D],Float64}, params[:order])
    reffeₚ = ReferenceFE(lagrangian, Float64, params[:order])


    V = TestFESpace(params[:model], reffeᵤ, conformity=:H1, dirichlet_tags=u_diri_tags)
    U = TransientTrialFESpace(V, u_diri_values)

    Q = TestFESpace(params[:model], reffeₚ, conformity=:H1, dirichlet_tags=p_diri_tags)
    P = TrialFESpace(Q, p_diri_values)

    Y = MultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])

    return V, U, P, Q, Y, X
end


"""
It creates the problem, here are implemented the VMS or SUPG system of equations
"""
function creation_op(params)
    if params[:method] == :VMS
        G, GG, gg = G_params(params[:Ω], params)
        if params[:linear]
            m, a, b = VMS_lin(G, GG, gg, params)
        else
            res, jac, jac_t = VMS(G, GG, gg, params)
        end
    elseif params[:method] == :SUPG
        h = h_param(params[:Ω], params[:D])
        if params[:linear]
            m, a, b = SUPG_lin(h, params)
        else
            res, jac, jac_t = SUPG(h, params)
        end
    else
        error("method $(params[:method]) not recognized\nuse :SUPG or :VMS")
    end
    if params[:linear]
        op = TransientAffineFEOperator(m, a, b, params[:X], params[:Y])
    else
        if params[:solver] == :petsc
            op = PETSC_TransientFEOperator(res, jac, jac_t, params[:X], params[:Y])
        elseif params[:solver] == :julia
            op = TransientFEOperator(res, jac, jac_t, params[:X], params[:Y])
        end
    end
    return op
end



"""
It creates the initial conditions. The default initial value for the pressure is 0.0 Pa. If the restart parameter is set to true, the initial conditons are set accordingly to the provided file
"""
function creation_initial_conditions(params::Dict{Symbol,Any})
    U0 = params[:U](0.0)
    P0 = params[:P](0.0)
    X0 = params[:X](0.0)
    if params[:restart] #for restarting from a file
        @info "restarting..."
        params[:u0] = restart_uh_field(params)

        p0 = restart_ph_field(params)

        if !(:p0 in keys(params))
            merge!(params, Dict(:p0 => p0))
        else
            params[:p0] = p0
        end

    end #end restart

    @info "interpolating uh"
    uh0 = interpolate_everywhere(params[:u0](0), U0)
    if params[:linear]
        ũ = uh0
        ũ_vector = create_ũ_vector(uh0.fields)
        merge!(params, Dict(:ũ => ũ, :ũ_vector => ũ_vector))
    end

    @info "interpolating ph"
    if :p0 in keys(params)              #p0 parameters for the TGV case, where there is an anyltic inital pressure
        ph0 = interpolate_everywhere(params[:p0](0), P0)
    else
        ph0 = interpolate_everywhere(0.0, P0)
    end

    xh0 = interpolate_everywhere([uh0, ph0], X0)
    if params[:restart] || params[:printinitial]
        writevtk(params[:Ω], "Initial_Condition", cellfields=["uh" => uh0, "ph" => ph0])
    end

    if params[:ode_method] == :AlphaMethod #Also the derivative is needed
        vu0 = (params[:D] == 2) ? VectorValue(0.0, 0.0) : VectorValue(0.0, 0.0, 0.0)
        vuh0 = interpolate_everywhere(vu0, U0)
        vph0 = interpolate_everywhere(0.0, P0)
        vxh0 = interpolate_everywhere([vuh0, vph0], X0)
        return (xh0, vxh0)

    elseif params[:ode_method] == :ThetaMethod
        return xh0
    end
end

"""
It creates parameters that are used for the ODE solution. 
In specific, if different ``\\theta`` are specified for velocity and pressure, are created vectors which relate each degree of freedom to the desired ``\\theta`` .
"""
function creation_ode_parameters(params::Dict{Symbol,Any})
    U0 = params[:U](0.0)
    P0 = params[:P](0.0)
    X0 = params[:X](0.0)
    @unpack θ, dt, D = params

    if typeof(θ) <: AbstractVector
        dtθ = dt .* θ

        if D == 2
            uhθt = interpolate_everywhere(VectorValue(θ[1], θ[1]), U0)
            uhθ = interpolate_everywhere(VectorValue(dtθ[1], dtθ[1]), U0)

        elseif D == 3
            uhθt = interpolate_everywhere(VectorValue(θ[1], θ[1], θ[1]), U0)
            uhθ = interpolate_everywhere(VectorValue(dtθ[1], dtθ[1], dtθ[1]), U0)
        end
        phθt = interpolate_everywhere(θ[2], P0)
        xhθt = interpolate_everywhere([uhθt, phθt], X0)

        θ_vec = get_free_dof_values(xhθt)

        phθ = interpolate_everywhere(dtθ[2], P0)
        xhθ = interpolate_everywhere([uhθ, phθ], X0)

        dtθ_vec = get_free_dof_values(xhθ)

        θ_params = [θ, θ_vec, dtθ, dtθ_vec]

    else
        θ_params = θ
    end
    return θ_params
    #merge!(params, Dict(:θ_params => θ_params))
end


"""
It creates the ODE solver using the parameters defined by the user
"""
function creation_ode_solver(params::Dict{Symbol,Any}) #In AlphaMethod case the xh0 is a tuple, first element is the initial condition, the second is the initial condition on the derivative
    θ_params = creation_ode_parameters(params)
    merge!(params, Dict(:θ_params => θ_params))
    @info "Ode parameters created"
    if params[:ode_method] == :AlphaMethod
        ode_solver = GeneralizedAlpha(params[:sys_solver], params[:dt], params[:ρ∞])

    elseif params[:ode_method] == :ThetaMethod
        ode_solver = ThetaMethod(params[:sys_solver], params[:dt], params[:θ_params])
        #ode_solver = ThetaMethod(params[:sys_solver], params[:dt], params[:θ_vec], params[:dtθ_vec])
    end

    sol_t = Algebra.solve(ode_solver, params[:op], params[:xh0], params[:t0], params[:tF])
    return sol_t, ode_solver
end

"""
It computes the solution, different for benchmark case or for getting the results.
"""
function compute_solution_(params::Dict{Symbol,Any})
    if params[:benchmark_mode]
        compute_solution_benchmark(params)
    else
        compute_solution(params)
    end
end

"""
It computes the solution, in the serial or parallel case.
"""
function compute_solution(params::Dict{Symbol,Any})
    @unpack backend = params
    simulation_filename = "$(params[:case])_$(params[:D])d"
    
        createpvd(params[:parts], simulation_filename) do pvd
            iterate_solution(params, pvd)
        end #end do
end

"""
It iterates the solution, here is where the solution is actually computed. At the end of each step the forces values are extracted if required
"""
function iterate_solution(params::Dict{Symbol,Any}, pvd)

    iter = 0
    for (xh_tn, tn) in params[:sol_t]
        iter = iter + 1
        println("iteration = $iter")
        uh_tn = xh_tn[1]
        ph_tn = xh_tn[2]
        ωh_tn = ∇ × uh_tn

        if !params[:benchmark_mode]
            pvd[tn] = createvtk(params[:Ω], joinpath("Results", "$(params[:case])_$(params[:D])d_$tn" * ".vtu"), cellfields=["uh" => uh_tn, "ph" => ph_tn, "wh" => ωh_tn])
        end

        write_forces(params, uh_tn, ph_tn, tn)

        update_linear!(params, uh_tn)
    end


end


"""
For benchmarking. It removes the compilation time for the first iteration and it does not print any results.
"""
function compute_solution_benchmark(params::Dict{Symbol,Any})
    @unpack tF, dt, sol_t, force_params, Ω, case = params
    @time (xh_tn, tn), state = Base.iterate(sol_t)
    uh_tn = xh_tn[1]
    update_linear!(params, uh_tn)
    @time while (tn < (tF - dt))
        (xh_tn, tn), state = Base.iterate(sol_t, state)
        uh_tn = xh_tn[1]
        update_linear!(params, uh_tn)
        println("Solution at time $tn")

    end #end while
end


"""
It creates the system solver which can be linear/non linear, julia/petsc using GAMG (only in petsc) or not.
"""
function create_system_solver(params::Dict{Symbol,Any})
    @unpack U, P = params
    sys_solver = :none

    if params[:solver] == :julia
        if params[:linear]
            sys_solver = LUSolver()
        else
            sys_solver = NLSolver(show_trace=params[:nls_trace], method=:newton, iterations=params[:nls_iter])
            @info "Using julia non linear solver"
        end

    elseif params[:solver] == :petsc
        if params[:linear]
            if params[:petsc_fieldsplit]
                FS = PETScFieldSplit([U(0), P(0)], ["vel", "pres"]; show_idx=false)
                sys_solver = PETScLinearSolver(mykspsetup, FS)
                @info "Using PETSc fieldsplit preconditioner"
            else
                sys_solver = PETScLinearSolver()
            end
            @info "Using PETSc linear solver"

        elseif params[:petsc_snes]
            if params[:petsc_fieldsplit]


                FS = PETScFieldSplit([U(0), P(0)], ["vel", "pres"]; show_idx=false)
                sys_solver = PETScNonlinearSolver(mysnessetup, FS)
                @info "Using PETSc fieldsplit preconditioner"

            else
                sys_solver = PETScNonlinearSolver()
            end
            @info "Using PETSc non linear solver"
        elseif params[:petsc_snes] == false
            ls = PETScLinearSolver()
            @info "Using PETSc linear solver"
            sys_solver = NLSolver(ls, show_trace=params[:nls_trace], method=:newton, iterations=params[:nls_iter])
            @info "Using Gridap non linear solver"
        end

    end
    return sys_solver
end

"""
Wrapper function that solve the case. Here the Non Linear Solver - nls is defined, in case the default julia sovler is adopted or PETSc
"""
function solve_case(params::Dict{Symbol,Any})


    sys_solver = create_system_solver(params)
    @info "system solver created"
    merge!(params, Dict(:sys_solver => sys_solver))

    xh0 = creation_initial_conditions(params)
    @info "initial conditions created"
    merge!(params, Dict(:xh0 => xh0))

    #xh0 = AccurateInitialCondition(ode_method, U, P, X, nls, xh0, op, t0, dt)
    #println("accurate initial conditions created")


    @time op = creation_op(params)
    @info "equations system created"
    merge!(params, Dict(:op => op))


    sol_t, ode_solver = creation_ode_solver(params)
    @info "ode solver created"
    merge!(params, Dict(:sol_t => sol_t, :ode_solver => ode_solver))

    #Print Jacobian and residual for the first time step, useful to verify that the matrices are created with the proper number of DoFs
    if params[:debug_mode] & false
        extract_matrix_vector(params)
    end

    params[:force_params] = initialize_force_params(params[:force_params], params[:case], params[:D], params[:ν])
    @info "force parameters created"

    """
    Profile.init()
    @profile @time compute_solution(params, params[:Ω])
    pprof()
    """
    @time compute_solution_(params)

    #@time compute_solution(params)
end
