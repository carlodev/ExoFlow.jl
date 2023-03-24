include("../../ExoFlow.jl")

using LinearAlgebra
using Statistics, Revise
params = Dict(
:u_in => 1.0,
:TI => 0.01,
:dt => 0.02,
:t0 => 0.0,
:Re_filename => joinpath("Channel", "Benchmark","Ret_395","Re_stress.xlsx" ), #"none", #
:dims => (:Y, ),
)

Vboxinfo = turbulence_box()

Re =  set_reynolds_stress(params)

sem_cache = SEMcache(Vboxinfo, params[:dt], params[:u_in], Vboxinfo.σ, Re, params[:t0], nothing)


Vboxinfo.N

x = VectorValue(1.0,0.5,0.0)

Nt = 0.1:0.1:10
U = zeros(length(Nt),3)
@time U_tmp = generation_u_fluct!(x, 0.01, sem_cache)

for i =1:1:length(Nt)
    t = Nt[i]
    U_tmp = generation_u_fluct!(x, t, sem_cache)
    U[i,1] = U_tmp[1]
    U[i,2] = U_tmp[2]
    U[i,3] = U_tmp[3]


end

