using ThermophysicalCalculation
using Base.Test
include ("test_IdealGasModel.jl")
include ("test_PengRobinsonModel.jl")
include ("test_Solver.jl")
forButane()
CompPolyVS_Hyper()
forAcetone()
forCarbonmonoxide()
