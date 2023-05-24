using ExoFlow
using LinearAlgebra
using Statistics, Revise, Gridap
using SyntheticEddyMethod
params = Dict(
:u_in => 1.0,
:TI => 0.01,
:dt => 0.01,
:t0 => 0.0,
:Re_filename => "none",
:dims => (:Y, ),
)

Vboxinfo = turbulence_box()

sem_cache.Eddies

Re =  set_reynolds_stress(params)

sem_cache = SEMcache(Vboxinfo, params[:dt], params[:u_in], Vboxinfo.σ, Re, params[:t0], nothing)

sem_cache.Re
Vboxinfo.N

x = VectorValue(0.0,0.0,0.0)

Nt = 0.0:params[:dt]:100
U = zeros(length(Nt),3)


for (i,t) in enumerate(Nt)
U_tmp = generation_u_fluct!(x, t, sem_cache)
U[i,1] = U_tmp[1]
U[i,2] = U_tmp[2]
U[i,3] = U_tmp[3]
println(i)
end
