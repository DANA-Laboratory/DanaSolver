using ThermophysicalCalculation
using Base.Test
using DanaTypes
#=
include ("test_Calls.jl")
include ("test_Solver.jl")
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
=#
include ("test_Analysis.jl")
