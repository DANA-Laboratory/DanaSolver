@test Solver.isconstantfactor([1.0,0.0,0.0,0.0,0.0,0.0])==true
function analysissystemofequations()
  println("**** analysis system of equations ****")
  exprs=[:(x+2*y+z^2),:(2*x-y*z),:(z-x+3*log(y)),:(x*y*z),:(x^y*log(z)),:(x+y),:(3*z-y+x),:(5.0+6.0)]
  re=Solver.analysis(exprs)
  println(exprs[6],"-->",re[1][1,:],'*',re[2],"=0")
  println(exprs[7],"-->",re[1][2,:],'*',re[2],"=0")
  @test length(re[3])==5 #index of noliexprs
  @test length(re[4])==8 #list of each expr args
  @test re[3][1]==1
  @test re[3][3]==3
  println ("list of args:")
  println (re[4])
end
