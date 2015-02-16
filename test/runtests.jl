using ThermophysicalCalculation
using Base.Test
using DanaTypes
#=
include ("test_Calls.jl")
include ("test_Analysis.jl")
allsyms()
analysissystemofequations()
include ("test_Solver.jl")
solvelinearidealgas()
testfindsystem()
solvenonlinearidealgas()
include ("test_IdealGasModel.jl")
forButane()
CompPolyVS_Hyper()
forwater()
forCarbonmonoxide()
=#

#include ("test_PengRobinsonDLModel.jl")
#testPRDLvriousunknowns()
include ("test_Analysis.jl")
#allsyms()
#analysissystemofequations()
testanalysis()
