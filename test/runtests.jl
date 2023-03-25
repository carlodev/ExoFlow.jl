using ExoFlow
using Test


include("driver.jl")
include("Turbulence_Settings.jl")

@testset "ExoFlow.jl" begin
   
    @test driver_test("TaylorGreen", :VMS, :ThetaMethod; linear = true )
    @test driver_test("TaylorGreen", :VMS, :ThetaMethod; linear = false)

    @test driver_test("TaylorGreen", :SUPG, :ThetaMethod; linear = true)
    @test driver_test("TaylorGreen", :SUPG, :ThetaMethod; linear = false )

    @test driver_test("TaylorGreen", :VMS, :AlphaMethod; linear = false)


    @test driver_test("LidDriven", :VMS, :ThetaMethod)
    @test driver_test("Cylinder", :VMS, :ThetaMethod; mesh_file = "Cylinder_2D.msh")

    @test driver_test("Airfoil", :VMS, :ThetaMethod )

end
