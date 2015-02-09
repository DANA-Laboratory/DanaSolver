using ThermophysicalCalculation
using Base.Test
using DanaTypes
include ("test_Solver.jl")
testfindsystem()
#=
include ("test_Calls.jl")
testforidealgasmodelwithcp()
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
include ("test_Analysis.jl")
=#
