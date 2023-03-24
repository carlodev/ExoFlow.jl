using Gridap
using Gridap.Algebra
using Gridap.Geometry
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Fields: meas

using GridapGmsh
using GridapDistributed
using GridapDistributed.CellData
#include("add_SEM_tag.jl")

model = GmshDiscreteModel("DU89_2D.msh")

#add_SEM_tag!(model)

L = 0.5
domain = (-L, L, -L, L)
partition = (params[:N], params[:N])
model = CartesianDiscreteModel(domain, partition, map=stretching)




writevtk(model, "NACA4412")


