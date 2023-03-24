

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


function physical_blasius(x, x_in, params)
@unpack u_in, D, ν, blasius_solution = params

f(η) = blasius_solution(η)[1]
fp(η) = blasius_solution(η)[2]
η(x,y) = y * (u_in / (ν*x))^0.5
u(x,y) =u_in * fp(η(x,y))
v(x,y) =0.5 * (ν* u_in/x)^0.5 *(η(x,y) * fp(η(x,y))- f(η(x,y)))

ub = u(x_in,x[2])
vb = v(x_in,x[2])
if D == 3
    return VectorValue(ub,vb,0.0)
else
    return VectorValue(ub,vb)
end

end


"""
using DifferentialEquations
using Gridap
using Parameters
blasius_solution = solve_blasius()
params = Dict(:balsius_solution => blasius_solution, :D =>2, :u_in => 1.0, :ν => 1e-5)
physical_blasius(Point(0.0,0.001), 0.005, params)
"""