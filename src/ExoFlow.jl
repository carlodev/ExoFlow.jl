"Package name"
module  ExoFlow

using Gridap
using Gridap.Algebra
using Gridap.Geometry
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Fields: meas
using Gridap.ODEs.ODETools
using Gridap.ODEs.TransientFETools
using Gridap.ODEs.ODETools: ThetaMethodNonlinearOperator
using Gridap.FESpaces: get_algebraic_operator


using GridapDistributed
using GridapDistributed.CellData

using GridapGmsh

using GridapPETSc
using GridapPETSc: PetscScalar, PetscInt, PETSC

using LineSearches: BackTracking, Static, MoreThuente
using FillArrays
using PartitionedArrays
using SparseArrays
using SparseMatricesCSR

# using DifferentialEquations
using JLD2
using MPI
using Random
using CSV
using DataFrames
using FileIO
# using BenchmarkTools
using Parameters
# using Profile
#using PProf
#using Revise
using Dates
using Revise
using XLSX
using Interpolations
using LinearAlgebra
using SyntheticEddyMethod

include("Main.jl")

#Equations
include(joinpath("Commons","Equations","Conservation.jl"))
include(joinpath("Commons","LinearUtilities.jl"))
include(joinpath("Commons","Equations","SUPGeq.jl"))
include(joinpath("Commons","Equations","VMSeq.jl"))
include(joinpath("Commons","PETSC_setups.jl"))

include(joinpath("Commons","Debug","Matrices_Vectors.jl"))
include(joinpath("Commons","Debug","Timing_tools.jl"))

include(joinpath("Commons","Tags","AddNewTags.jl"))

include(joinpath("Commons","SEM","SEM_interface.jl"))
include(joinpath("Commons","SEM","Turbulence_Utilities.jl"))


include(joinpath("Commons","ComputeForces.jl"))
include(joinpath("Commons","CommonProcedures.jl"))
include(joinpath("Commons","Errors_Warnings.jl"))
include(joinpath("Commons","Init_params.jl"))
include(joinpath("Commons","Restart.jl"))

#Cases
include(joinpath("Channel","Channel.jl"))
include(joinpath("LidDriven","LidDriven.jl"))
include(joinpath("Cylinder","Cylinder.jl"))
include(joinpath("Airfoil","Airfoil.jl"))
include(joinpath("Plate","Plate.jl"))
include(joinpath("Plate","Blasius.jl"))
include(joinpath("TaylorGreen","TaylorGreen.jl"))

export hf_gen!
export add_centre_tag!
export printmodel
export creation_fe_spaces
export creation_op
export creation_initial_conditions
export creation_ode_solver
export compute_solution
export iterate_solution
export solve_case

export petsc_options
export initialize_parameters
export set_reynolds_stress
end # module