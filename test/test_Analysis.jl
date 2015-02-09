@test Solver.isconstantfactor([1.0,0.0,0.0,0.0,0.0,0.0])==true
exprs=[:(x+2*y+z^2),:(2*x-y*z),:(z-x+3*log(y)),:(x*y*z),:(x^y*log(z)),:(x+y),:(3*z-y+x),:(5.0+6.0)]
re=Solver.analysis(exprs)
println(exprs[6],"-->",re[1][1,:],'*',re[2],"=0")
println(exprs[7],"-->",re[1][2,:],'*',re[2],"=0")
@test length(re[3])==5
@test length(re[4])==5
@test re[3][1](0.0,0.0,2.0)==4
@test re[3][3](1.0,1.0,1.0)==0.0
re=Solver.analysis(exprs,true)
@test length(re[3])==8
@test length(re[4])==8
@test re[3][1](0.0,0.0,2.0)==4
@test re[3][3](1.0,1.0,1.0)==0.0
@test re[3][6](1.0,1.0)==2.0
