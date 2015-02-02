function testforidealgasmodelwithcp()
  DNIdel=DANAIdealGasEos()
  DNIdel.P=12.0
  DNIdel.T=120.0
  DNIdel.CASNO="95-63-6" #1,2,4-Trimethylbenzene
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly","1,2,4-Trimethylbenzene")
  setEquationFlow(DNIdel)
  a=Solver.replace(DNIdel)
  println(a[5])
  a=map(Calculus.simplify,a)
  vals,vars=Solver.analysis(a)
  println(a)
  println(vals)
  println(vars)
  rVls=Solver.rref(vals)
  DNIdel.v=-1*last(rVls[1,:])
  DNIdel.Cp=-1*last(rVls[2,:])
  println(rVls)
  println(vars)
end
function testupdate()
  DNIdel=DANAIdealGasEos()
  DNIdel.P=12.0
  DNIdel.T=120.0
  DNIdel.CASNO="95-63-6" #1,2,4-Trimethylbenzene
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly","1,2,4-Trimethylbenzene")
  setEquationFlow(DNIdel)
  rVls,vars=Solver.solvelinear(DNIdel)
  Solver.update!(DNIdel,rVls,vars)
  a=Solver.replace(DNIdel)
  println(map(eval,a))
end
function testidealgas()
  ######Temprature is undef#######
  DNIdel=DANAIdealGasEos()
  DNIdel.P=2000.0
  DNIdel.v=4000.0
  DNIdel.CASNO="95-63-6" #1,2,4-Trimethylbenzene
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly","1,2,4-Trimethylbenzene")
  setEquationFlow(DNIdel)
  somethingUpdated,fullDetermined=Solver.slstsubnfd!(DNIdel)
  dump(DNIdel)
  a=Solver.replace(DNIdel)
  println(a)
end
function testidealgasfornonlinearsolver()
  ######Temprature is undef#######
  DNIdel=DANAIdealGasEos()
  DNIdel.P=2000.0
  DNIdel.Cp=629657.0
  DNIdel.CASNO="95-63-6"
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly","1,2,4-Trimethylbenzene")
  setEquationFlow(DNIdel)
  somethingUpdated=true
  fullDetermined=false
  while (somethingUpdated && !fullDetermined)
    rVls,vars=Solver.solvelinear(DNIdel)
    println("************one solution done************")
    somethingUpdated,fullDetermined=Solver.update!(DNIdel,rVls,vars)
  end
  dump(DNIdel)
end
function testidealgasfornonlinearsolver2()
  ######Temprature is undef###using nolinear solver for only an equation####
  DNIdel=DANAIdealGasEos()
  DNIdel.P=2000.0
  DNIdel.Cp=629657.0
  DNIdel.CASNO="95-63-6"
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly","1,2,4-Trimethylbenzene")
  setEquationFlow(DNIdel)
  somethingUpdated=true
  fullDetermined=false
  nonliVars::Array{Set{String},1}=Array(Set{String},0)
  nonliFuns::Array{Function,1}=Array(Function,0)
  while (somethingUpdated && !fullDetermined)
    while (somethingUpdated && !fullDetermined)
      rVls,vars,nonliFuns,nonliVars=Solver.solvelinear(DNIdel)
      println("************one linear solution done************")
      somethingUpdated,fullDetermined=Solver.update!(DNIdel,rVls,vars)
    end
    println("fullDetermined=",fullDetermined)
    if !fullDetermined
      i=1
      fullDetermined=true
      while (i<=length(nonliFuns))
        if length(nonliVars[i])==1
          result=Roots.fzero(nonliFuns[i],[0,typemax(Int64)])
          Solver.setfield!(DNIdel,[nonliVars[i]...][1],result)
          println("************one nonlinear solution done************")
          somethingUpdated=true
        else
          fullDetermined=false
        end 
        i=i+1
      end
    end
    println("fullDetermined=",fullDetermined)
  end
  println("somethingUpdated=",somethingUpdated)
  println("fullDetermined=",fullDetermined)
  dump(DNIdel)
end
function testsolve()
  DNIdel=DANAIdealGasEos()
  DNIdel.P=2000.0
  DNIdel.Cp=629657.0
  DNIdel.CASNO="95-63-6"
  DNIdel.usePolynomialEstimationOfCp=true
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpPoly","1,2,4-Trimethylbenzene")
  setEquationFlow(DNIdel)
  somethingUpdated,fullDetermined,noliTrys,nonlTrys=Solver.sliponl!(DNIdel)
  println("somethingUpdate=",somethingUpdated," fullDetermined=",fullDetermined," noliTrys=",noliTrys," nonlTrys=",nonlTrys)
end
