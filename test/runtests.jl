using ExoFlow
using Test

@testset "ExoFlow.jl" begin
    include("driver.jl")
    @test driver_test("TaylorGreen", :VMS, :ThetaMethod; linear = true )
    @test driver_test("TaylorGreen", :VMS, :ThetaMethod; linear = false)

    @test driver_test("TaylorGreen", :SUPG, :ThetaMethod; linear = true)
    @test driver_test("TaylorGreen", :SUPG, :ThetaMethod; linear = false )

    @test driver_test("TaylorGreen", :VMS, :AlphaMethod; linear = false)


    @test driver_test("LidDriven", :VMS, :ThetaMethod )
    @test driver_test("Airfoil", :VMS, :ThetaMethod )

end
