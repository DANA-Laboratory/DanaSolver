function solvelinearexprforvar()
  println("**** test solve linear expression for var with fac ****")
  @test Solver.solveexpr!(:(5x+y),:x,5.0) == :(-1y/5)
  @test Solver.solveexpr!(:(sin(z)*y+5x+y^2+z),:x,5.0) == :(-1(sin(z)*y+y^2+z)/5)
end
function solvelinearidealgas()
  println("**** solve linear ideal gas ****")
  println("****** 1-step by step ******")
  DNIdel=DANAIdealGasEos()
  DNIdel.P=12.0
  DNIdel.T=120.0
  DNIdel.name="1,2,4-Trimethylbenzene" #95-63-6 
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly",DNIdel.name)
  setEquationFlow(DNIdel)
  a=Solver.replace(DNIdel)
  a=map(Calculus.simplify,a)
  vals,vars=Solver.analysis(a)
  rVls=Solver.rref(vals)
  @test_approx_eq 83144.621 -1*last(rVls[1,:])
  @test_approx_eq 78910.7999999 -1*last(rVls[2,:])
  println("****** 2-using linear solver ******")
  DNIdel=DANAIdealGasEos()
  DNIdel.P=12.0
  DNIdel.T=120.0
  DNIdel.name="1,2,4-Trimethylbenzene" #95-63-6 
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly",DNIdel.name)
  setEquationFlow(DNIdel)
  rVls,vars=Solver.solvelinear(DNIdel)
  Solver.update!(DNIdel,rVls,vars)
  a=Solver.replace(DNIdel)
  @test zeros(8) == (map(eval,a))
  println("****** 3-using slstsubnfd! ******")
  DNIdel=DANAIdealGasEos()
  DNIdel.P=12.0
  DNIdel.T=120.0
  DNIdel.name="1,2,4-Trimethylbenzene" #95-63-6 
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly",DNIdel.name)
  somethingUpdated,fullDetermined,nonliExprIndx,args,equations,noTrys=Solver.slstsubnfd!(DNIdel)
  @test somethingUpdated == true
  @test fullDetermined == true
  @test nonliExprIndx == []
  @test length(args) == 8 
  @test length(equations) == 8 
  println("solution done by $noTrys trys")
  @test_approx_eq 83144.621 DNIdel.v
  @test_approx_eq 78910.7999999 DNIdel.Cp
  println("****** 4-using solve ******")
  DNIdel=DANAIdealGasEos()
  DNIdel.P=12.0
  DNIdel.T=120.0
  DNIdel.name="1,2,4-Trimethylbenzene" #95-63-6 
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly",DNIdel.name)
  somethingUpdated,fullDetermined,noliTrys,nonlTrys=solve!(DNIdel)
  @test somethingUpdated == true
  @test fullDetermined == true
  @test noliTrys == 1
  @test nonlTrys == 0
  @test_approx_eq 83144.621 DNIdel.v
  @test_approx_eq 78910.7999999 DNIdel.Cp
end
function solvenonlinearidealgas()
  println("**** solve nonlinear ideal gas ****")
  println("****** 1-using fzero ******")
  ######Temprature is undef#######
  DNIdel=DANAIdealGasEos()
  DNIdel.P=2000.0
  DNIdel.Cp=629657.0
  DNIdel.name="1,2,4-Trimethylbenzene" #"95-63-6"
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly",DNIdel.name)
  somethingUpdated,fullDetermined,nonliExprIndx,args,equations,noTrys=Solver.slstsubnfd!(DNIdel)
  @test somethingUpdated == false
  @test fullDetermined == false
  println("linear solution done by $noTrys trys")
  println("model has ",length(nonliExprIndx)," nolinear equations->")
  println(getindex(equations,nonliExprIndx))
  varIndex,allVars,eqIndex=Solver.findsystem(args)
  println(varIndex,allVars,eqIndex)
  if length(eqIndex)==1
    fun=Solver.exprTofunction(equations[eqIndex[1]],allVars[varIndex[1]])
    result=Roots.fzero(fun,[0,typemax(Int64)])
    @test_approx_eq result 962.1783117595058 
    Solver.setfield!(DNIdel,allVars[varIndex[1]][1],result)
  end
  @test_approx_eq DNIdel.T 962.1783117595058
  println("****** 2-using sliponl! ******")
  DNIdel=DANAIdealGasEos()
  DNIdel.P=2000.0
  DNIdel.Cp=629657.0
  DNIdel.name="1,2,4-Trimethylbenzene" #"95-63-6"
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly",DNIdel.name)
  somethingUpdated,fullDetermined,noliTrys,nonlTrys=Solver.sliponl!(DNIdel)
  @test somethingUpdated == true
  @test fullDetermined == true
  @test_approx_eq DNIdel.T 962.1783117595058
  @test_approx_eq DNIdel.ICpDT 2.7479202136168545e8 
  @test_approx_eq DNIdel.ICpOnTDT 697723.1775663692 
  println ("solution done using $noliTrys linear plus $nonlTrys nonlinear trys")
  println("****** 3-using solve! ******")
  DNIdel=DANAIdealGasEos()
  DNIdel.P=2000.0
  DNIdel.Cp=629657.0
  DNIdel.name="1,2,4-Trimethylbenzene" #"95-63-6"
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly",DNIdel.name)
  somethingUpdated,fullDetermined,noliTrys,nonlTrys=solve!(DNIdel)
  @test somethingUpdated == true
  @test fullDetermined == true
  @test_approx_eq DNIdel.T 962.1783117595058
  @test_approx_eq DNIdel.ICpDT 2.7479202136168545e8 
  @test_approx_eq DNIdel.ICpOnTDT 697723.1775663692 
  println ("solution done using $noliTrys linear plus $nonlTrys nonlinear trys")
end
function testfindsystem()
  println("**** find system to solve ****")
  args=[Set(String["c","d","c"]),Set(String["b","c","i"]),Set(String["d","g","h"]),Set(String["d","c","i"]),Set(String["a","d","k"]),Set(String["b","a","c"]),Set(String["a","b","c"]),Set(String["a","d"]),Set(String["d","e"]),Set(String["e","f"]),Set(String["f","b","e"]),Set(String["e","k"])]
  println("****** test for a system of ",length(args)," equations-> ******")
  println(args)
  varindex,allvars,eqindex=Solver.findsystem(args)
  i=1
  println("equation number    |    in tems of ")
  for ai in varindex
    println("       ",eqindex[i],"           | ",allvars[ai])
    i+=1
  end
end
