

"""
    solve_blasius()

It solves Blasius equations and provides 3 values: 
f
``\\dfrac{\\partial f}{\\partial \\eta}``
``\\dfrac{\\partial^2  f}{\\partial \\eta^2}``  
# Example
```julia-repl
julia> blasius_solution = solve_blasius()
julia> blasius_solution(0.5)
3-element Vector{Float64}:
 0.04149305600494741
 0.16588619765183715
 0.3309128332084983
 ```
"""
function solve_blasius()

    function balsius_fun!(du, u, p, η)
    f, f1, f2 = u
    du[1] = f1
    du[2] = f2
    du[3] = -0.5*f2*f
end

    function bc_blasius!(residual, u, p, η) # u[1] is the beginning of the time span, and u[end] is the ending
        residual[1] = u[1][1]   # the solution at the beginning of the time span should be -pi/2
        residual[2] = u[1][2]  # the solution at the end of the time span should be pi/2
        residual[3] = u[end][2] - 1.0
    end
        ηspan = (0,8)

    blasius_problem = TwoPointBVProblem(balsius_fun!, bc_blasius!, [0.0, 0.0, 0.3321], ηspan)
    balsius_solution = DifferentialEquations.solve(blasius_problem, Shooting(Vern7()), dt=0.001) #dt in reality dη
    return balsius_solution
end

"""
    physical_blasius(x::Point, x_in::Float64, params::Dict{Symbol,Any})

It computes (u,v) or (u,u,w), the velocity in 2 or 3 direction, using the Blasius' equations.
It computes the boundary layer properties having x_in has development length.
```julia
using DifferentialEquations
using Gridap
using Parameters
blasius_solution = solve_blasius()
params = Dict(:blasius_solution => blasius_solution, :D =>2, :u_in => 1.0, :ν => 1e-5)

physical_blasius(Point(0.01, 0.01), 1.0, params)
```

"""
function physical_blasius(x::Point, x_in::Float64, params::Dict{Symbol,Any})
@unpack u_in, D, ν = params

blasius_solution = solve_blasius()
function f(η)
    if η <8
    return blasius_solution(η)[1]
    else
        β = 1.720787657
        return η - β
    end
end

function fp(η)
    if η <8
        return blasius_solution(η)[2]
        else
            return 1
        end
end

η(x,y) = y * (u_in / (ν*x))^0.5
u(x,y) = u_in * fp(η(x,y))
v(x,y) =0.5 * (ν* u_in/x)^0.5 *(η(x,y) * fp(η(x,y))- f(η(x,y)))

ub = u(x_in,x[2])
vb = v(x_in,x[2])
if D == 3
    return VectorValue(ub,vb,0.0)
else
    return VectorValue(ub,vb)
end

end
function physical_blasius(x::Point, params::Dict{Symbol,Any})
x_in = x[1]
physical_blasius(x, x_in, params)
end
