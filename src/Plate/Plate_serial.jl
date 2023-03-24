include("../ExoFlow.jl")
include("BoundaryConditions.jl")

function stretching(x::Point)
    m = zeros(length(x))
    m[1] = x[1]
    m[2] = x[2].^2
    Point(m)
  end
function petsc_gamg_options()
    "-ksp_type gmres -pc_type gamg -pc_gamg_type agg -ksp_rtol 1.0e-9"
end
"""
-ksp_type gmres -ksp_rtol 1.0e-06 -ksp_atol 0.0 -pc_type gamg -pc_gamg_type agg \
      -mg_levels_esteig_ksp_type gmres -mg_coarse_sub_pc_type lu \
      -mg_coarse_sub_pc_factor_mat_ordering_type nd -pc_gamg_process_eq_limit 50 \
      -pc_gamg_square_graph 9 pc_gamg_agg_nsmooths 1
      """
parameters = Dict(
    :N => 32,
    :D => 2, #Dimension
    :order => 1, 
    :t0 =>0.0,
    :dt => 0.01, #0.01
    :tF => 10.0,
    :t_endramp => 5.0,

    :case => "Plate",
    :solver => :julia,
    :method => :VMS,
    :ode_method => :AlphaMethod, 
    :ρ∞ => 0.8, #0.01
    :Re => 500_000,
    :c => 1.0, #chord lenght [m], used for naca (=1), cylinder (=0.1), liddriven = (1), 0.5
    :u_in => 10.0,  # =1.0 for lid driven 
    :periodic => false,
    :printmodel => true,

    :mesh_file => "NACA0012_3D_improved.msh",
    :Cᵢ => [10, 1], 
    :options => petsc_gamg_options(), 
    :ksp_monitor => false, #do not print every iteration
    :nls_trace =>true,

    :ν => 1e-5, #channel = 0.0001472, 
    :ρ => 1.0, #kg/m3 density
    :body_force => 0.0, #channel = 0.00337204

    :np_x => 8, #number of processors in X
    :np_y => 8, #number of processors in Y
    :np_z => 2, #number of processors in Z

    :restart => false,
    :restart_file => "Airfoil/DU89_319.csv"
)

params = initialize_parameters(parameters)
params[:ν] = 1e-5 #m2/s 

Lx = 3
Ly = 1
Nx = 250
Ny = 400
    

domain = (0, Lx, 0, Ly)
partition = (Nx, Ny)
model = CartesianDiscreteModel(domain, partition,map=stretching)


params[:printmodel]
printmodel(params, model)
hf_gen!(params)


u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u0 = bc_plate(params) 
merge!(params, Dict(:u0 => u0, :model => model))

U,P,Y,X = creation_fe_spaces(params, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)

degree = 4*params[:order]
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

new_dict = Dict(:U => U,
    :P => P,
    :X => X,
    :Y => Y,
    :Ω => Ω,
    :dΩ => dΩ,
    :degree => degree,
    :force_params => nothing)
merge!(params, new_dict)


solve_case(params)
