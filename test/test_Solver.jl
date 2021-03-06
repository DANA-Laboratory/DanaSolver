function testwhatif()
  @test Solver.whatifcombineeqs(Set{String}(),Set{String}(),Set{String}(),Set{String}()) == []
  @test Solver.whatifcombineeqs(Set{String}(["x","y"]),Set{String}(["z"]),Set{String}(["y"]),Set{String}(["x"])) == [Set{String}(),Set{String}(["z","x"]),Set{String}(),Set{String}(["z","y"])]
  @test Solver.whatifcombineeqs(Set{String}(["x","y"]),Set{String}(["z"]),Set{String}(["x","y"]),Set{String}(["z"])) == [Set{String}(["y"]),Set{String}(["z"]),Set{String}(["x"]),Set{String}(["z"])]
  @test Solver.whatifcombineeqs(Set{String}(["x","y"]),Set{String}(["z"]),Set{String}(["z"]),Set{String}(["x"])) == [Set{String}(),Set{String}(["z","y"]),Set{String}(["y"]),Set{String}(["x"])]
end
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
  eqs=Solver.Equations(map(Solver.Equation,a))
  Solver.analysis(eqs)
  vals,vars=eqs.facts,eqs.terms
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
  rVls,eqs=Solver.solvelinear(DNIdel)
  Solver.update!(DNIdel,rVls,eqs.terms)
  a=Solver.replace(DNIdel)
  @test zeros(8) == (map(eval,a))
  println("****** 3-using slstsubnfd! ******")
  DNIdel=DANAIdealGasEos()
  DNIdel.P=12.0
  DNIdel.T=120.0
  DNIdel.name="1,2,4-Trimethylbenzene" #95-63-6
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly",DNIdel.name)
  somethingUpdated,fullDetermined,eqs,noTrys=Solver.slstsubnfd!(DNIdel)
  @test somethingUpdated == true
  @test fullDetermined == true
  @test eqs.indexnonliexs == []
  @test length(eqs.terms) == 9
  @test length(eqs.exs) == 8
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
  somethingUpdated,fullDetermined,eqs,noTrys=Solver.slstsubnfd!(DNIdel)
  @test somethingUpdated == false
  @test fullDetermined == false
  println("linear solution done by $noTrys trys")
  println("model has ",length(eqs.indexnonliexs)," nolinear equations->")
  println(getindex(eqs.exs,eqs.indexnonliexs))
  varIndex,allVars,eqIndex=Solver.findsystem(eqs)
  println(varIndex,allVars,eqIndex)
  if length(eqIndex)==1
    fun=Solver.exprTofunction(eqs.exs[eqIndex[1]].ex,allVars[varIndex[1]])
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
  println("solution done using $noliTrys linear plus $nonlTrys nonlinear trys")
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
  println("solution done using $noliTrys linear plus $nonlTrys nonlinear trys")
end
function testfindsystem()
  println("**** find system to solve ****")
  args=[Set(String["c","d","c"]),Set(String["b","c","i"]),Set(String["d","g","h"]),Set(String["d","c","i"]),Set(String["a","d","k"]),Set(String["b","a","c"]),Set(String["a","b","c"]),Set(String["a","d"]),Set(String["d","e"]),Set(String["e","f"]),Set(String["f","b","e"]),Set(String["e","k"])]
  argsnli=[Set(String["c"]),Set(String["b"]),Set(String["g","h"]),Set(String["d","i"]),Set(String[]),Set(String["a","c"]),Set(String[]),Set(String["d"]),Set(String["e"]),Set(String["e"]),Set(String[]),Set(String["e"])]
  println("****** test for a system of ",length(args)," equations-> ******")
  println(args)
  eqs=Array(Solver.Equation,12)
  i=1;
  while i<=length(args) ; eqs[i]=Solver.Equation(:()); eqs[i].termall=args[i]; eqs[i].termnonli=argsnli[i]; i+=1; end;
  j=1
  varindex,allvars,eqindex=false,false,false
  while varindex==false
    varindex,allvars,eqindex=Solver.findsystem(Solver.Equations(eqs),j)
    j+=1
  end
  i=1
  println("number of nolinear terms=$j")
  println("equation number    |    in tems of ")
  for ai in varindex
    println("       ",eqindex[i],"           | ",allvars[ai])
    i+=1
  end
end
