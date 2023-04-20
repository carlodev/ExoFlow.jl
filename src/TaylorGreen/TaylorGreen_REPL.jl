using Revise
using ExoFlow
using Gridap, GridapPETSc
using Profile
using ProfileView
using PProf
include("SpaceConditions.jl")
include("AnalyticalSolution.jl")
N = 128
D = 2
order = 1
t0 = 0.0
dt = 0.01
tF = 0.10
diameter = 0.5
θ = [1.0,1.0]
Vs = 1.0 #1[m/s]swirling speed
Ua = 0.3 #0.3 [m/s]convective velocity in x
Va = 0.2 #0.2 [m/s]convective velocity in y
ν = 0.001 #0.001 m2/s 
body_force = 0.0
Cᵢ=[4,36]

options ="-pc_type asm -sub_pc_type ilu -ksp_type gmres -ksp_gmres_restart 30  -snes_converged_reason -ksp_converged_reason -ksp_error_if_not_converged true"
#options ="-pc_type none -ksp_type gmres -ksp_gmres_restart 30  -snes_converged_reason -ksp_converged_reason -ksp_error_if_not_converged true"

# options = options * " -log_view" * " -ksp_monitor_short"

domain = (-diameter, diameter, -diameter, diameter)
partition = (N,N)
model = CartesianDiscreteModel(domain, partition; isperiodic=(true, true))
velocity, pa, ωa = analytical_solution(diameter, Vs, Ua, Va, ν)
params = Dict(:D=> D, :order => order, :body_force=> body_force, :θ=>θ, :linear=>true, :method=>:VMS)
ExoFlow.hf_gen!(params)
hf = params[:hf]
V, Q, U, P, Y, X, model = CreateTGSpaces(model, params, pa) #We update model with the new label of the center point

degree = 4 * order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)
uh0 = interpolate_everywhere(velocity(0), U(0))
ph0 = interpolate_everywhere(pa(0), P(0))
xh0 = interpolate_everywhere([uh0,ph0], X(0))
ũ_vector = create_ũ_vector(uh0)
ũ = uh0

Cᵢ1 = Cᵢ[1]
Cᵢ2 = Cᵢ[2]
merge!(params, Dict(:Ω=>Ω, :dΩ=>dΩ, :debug_mode=>false, :ν=>ν, :dt=>dt, :ũ => ũ, :ũ_vector => ũ_vector, :Cᵢ=>Cᵢ, :steady=>false, :X => X, :Y => Y))
merge!(params, Dict(:U=> U, :P => P))
ũ_coeff = [2.1875, -2.1875, 1.3125, -0.3125]
merge!(params, Dict(:ũ_coeff => ũ_coeff))

#Stabilization
G, GG, gg = G_params(Ω, Dict(:D=>2, :debug_mode=>false))
G
uh0.free_values

function τm(uu,G, GG)
  τ₁ = Cᵢ1 * (2 / dt)^2 #Here, you can increse the 2 if CFL high
  τ₃ = Cᵢ2 * (ν^2 * GG)

  val(x) = x
  function val(x::Gridap.Fields.ForwardDiff.Dual)
    x.value
  end
  D = length(uu)
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

  if iszero(norm(uu_new))
    return (τ₁ .+ τ₃) .^ (-1 / 2)
  end

  τ₂ = uu_new ⋅ G ⋅ uu_new
  return (τ₁ .+ τ₂ .+ τ₃) .^ (-1 / 2)
end

#Stabilization paramameter continuity
function τc(uu,G, GG,gg)
  return 1 / (τm(uu, G, GG) ⋅ gg)
end

#Equations

θvp = 1
aRm, aRc = conservations_equations(params)
aTRm(u, p) = τm ∘ (ũ,G, GG) * aRm(u, p)
bTRm(u, p) = (τm ∘ (ũ,G, GG)) * u
bbTRm(u, p) = (τm ∘ (ũ,G, GG)) * ũ

fTRm(t) = τm ∘ (ũ, G, GG) * hf(t)


m(t, (u, p), (v, q)) =  ∫(u ⋅ v)dΩ + 
∫((ũ ⋅ ∇(v) + (θvp)*∇(q)) ⊙ bTRm(u, p))dΩ + 
 ∫((ũ ⋅ (∇(v))') ⊙ bTRm(u, p))dΩ #+
#  -0.5 * ∫((∇(v)) ⊙ (outer(bTRm(u, p), bbTRm(u, p))))dΩ +
#  -0.5 * ∫((∇(v)) ⊙ (outer(bbTRm(u, p), bTRm(u, p))))dΩ


aᴳ(t, (u, p), (v, q)) = ∫((ũ ⋅ ∇(u)) ⋅ v)dΩ - ∫((∇ ⋅ v) * p)dΩ + ∫((q * (∇ ⋅ u)))dΩ + ν * ∫(∇(v) ⊙ ∇(u))dΩ 

a_SUPG(t, (u, p), (v, q)) = ∫((ũ ⋅ ∇(v) + ∇(q)) ⊙ aTRm(u, p))dΩ + ∫((∇ ⋅ v) ⊙ (τc ∘ (ũ,G, GG,gg) * aRc(u)))dΩ
a_VMS1(t,(u, p), (v, q)) = ∫((ũ ⋅ (∇(v))') ⊙ aTRm(u, p))dΩ
a_VMS2(t, (u, p), (v, q)) = -0.5 * ∫((∇(v)) ⊙ (outer(aTRm(u, p), bbTRm(u, p))))dΩ -0.5* ∫((∇(v)) ⊙ (outer(bbTRm(u, p), aTRm(u, p))))dΩ

a(t, (u, p), (v, q)) = aᴳ(t, (u, p), (v, q)) +  a_SUPG(t, (u, p), (v, q)) + a_VMS1(t,(u, p), (v, q)) #+ a_VMS2(t,(u, p), (v, q))
b(t, (v, q)) =∫(hf(t) ⋅ v)dΩ + ∫((ũ ⋅ ∇(v) + ∇(q)) ⋅  fTRm(t))dΩ + ∫((ũ ⋅ (∇(v))') ⋅ fTRm(t))dΩ
op = TransientAffineFEOperator(m, a, b,X,Y)

θ_params = creation_ode_parameters(params)

using Gridap
tF = 5*dt
GridapPETSc.with(args=split(options)) do
  ls = PETScLinearSolver()
  # ls = LUSolver()
  ode_solver = ThetaMethod(ls, dt, θ_params)
  sol_t = Gridap.Algebra.solve(ode_solver, op, xh0, t0, tF)
  
  (xh_tn, tn), state = Base.iterate(sol_t)

    uh_tn = xh_tn[1]
    update_linear!(params, uh_tn)
    # ProfileView.@profview # Profile.Allocs.@profile sample_rate = 0.001 
   @time while (tn < (tF - dt))
       (xh_tn, tn), state = Base.iterate(sol_t, state)
        uh_tn = xh_tn[1]
         update_linear!(params, uh_tn)
        println("Solution at time $tn")

    end #end while

 end #GridapPetscdo



# PProf.Allocs.pprof(from_c = false)