using ThermophysicalCalculation
using Base.Test
using DanaTypes
include ("test_Analysis.jl")
analysissystemofequations()
include ("test_Solver.jl")
solvelinearidealgas()
testfindsystem()
solvenonlinearidealgas()
#=
include ("test_Calls.jl")
testupdate()
testidealgas()
testidealgasfornonlinearsolver()
testidealgasfornonlinearsolver2()
testsolve()
include ("test_IdealGasModel.jl")
forButane()
CompPolyVS_Hyper()
forAcetone()
forCarbonmonoxide()
include ("test_PengRobinsonDLModel.jl")
testMoreThanOneNonLinear()
=#
