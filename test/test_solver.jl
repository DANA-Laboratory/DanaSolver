reload ("HelperEquation.jl")
reload ("IdealGasEos.jl")
reload ("Tables.jl")
module test
  using HelperEquation
  using IdealGasEos
  using Tables
  using Roots
  function testforIdealGasModelWithCp()
    DNIdel=DANAIdealGasEos()
    DNIdel.P=12.0
    DNIdel.T=120.0
    DNIdel.CASNO="95-63-6"
    DNIdel.usePolynomialEstimationOfCp=true
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getValueForCasNo("C0Poly","95-63-6")
    setEquationFlow(DNIdel)
    a=replace(DNIdel)
    println(a[5])
    a=map(Calculus.simplify,a)
    vals,vars=analysis(a)
    println(a)
    println(vals)
    println(vars)
    rVls=rref(vals)
    DNIdel.v=-1*last(rVls[1,:])
    DNIdel.Cp=-1*last(rVls[2,:])
    println(rVls)
    println(vars)
  end
  function testUpdate()
    DNIdel=DANAIdealGasEos()
    DNIdel.P=12.0
    DNIdel.T=120.0
    DNIdel.CASNO="95-63-6"
    DNIdel.usePolynomialEstimationOfCp=true
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getValueForCasNo("C0Poly","95-63-6")
    setEquationFlow(DNIdel)
    rVls,vars=solve(DNIdel)
    update!(DNIdel,rVls,vars)
    a=replace(DNIdel)
    println(map(eval,a))
  end
  function testIDealGas()
    ######Temprature is undef#######
    DNIdel=DANAIdealGasEos()
    DNIdel.P=2000.0
    DNIdel.v=4000.0
    DNIdel.CASNO="95-63-6"
    DNIdel.usePolynomialEstimationOfCp=true
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getValueForCasNo("C0Poly","95-63-6")
    setEquationFlow(DNIdel)
    somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      rVls,vars=solve(DNIdel)
      somthingUpdated,fullDetermined=update!(DNIdel,rVls,vars)
    end
    dump(DNIdel)
    #a=replace(DNIdel)
    #println(a)
  end
  function testIDealGasForNonlinearSolver()
    ######Temprature is undef#######
    DNIdel=DANAIdealGasEos()
    DNIdel.P=2000.0
    DNIdel.Cp=629657.0
    DNIdel.CASNO="95-63-6"
    DNIdel.usePolynomialEstimationOfCp=true
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getValueForCasNo("C0Poly","95-63-6")
    setEquationFlow(DNIdel)
    somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      rVls,vars=solve(DNIdel)
      println("************one solution done************")
      somthingUpdated,fullDetermined=update!(DNIdel,rVls,vars)
    end
    dump(DNIdel)
  end
  function testIDealGasForNonlinearSolver2()
    ######Temprature is undef###using nolinear solver for only an equation####
    DNIdel=DANAIdealGasEos()
    DNIdel.P=2000.0
    DNIdel.Cp=629657.0
    DNIdel.CASNO="95-63-6"
    DNIdel.usePolynomialEstimationOfCp=true
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getValueForCasNo("C0Poly","95-63-6")
    setEquationFlow(DNIdel)
    somthingUpdated=true
    fullDetermined=false
    nonliFuns::Array{Function,1}=Array(Function,0)
    nonliVars::Array{Array{String,1},1}=Array(Array{String,1},0)
    while (somthingUpdated && !fullDetermined)
      while (somthingUpdated && !fullDetermined)
        rVls,vars,nonliFuns,nonliVars=solve(DNIdel)
        println("************one linear solution done************")
        somthingUpdated,fullDetermined=update!(DNIdel,rVls,vars)
      end
			println(fullDetermined)
      if !fullDetermined
        i=1
        fullDetermined=true
        while (i<=length(nonliFuns))
          if length(nonliVars[i])==1
            HelperEquation.setfield(DNIdel,nonliVars[i][1],Roots.fzero(nonliFuns[i],[0,typemax(Int64)]))
            println("************one nonlinear solution done************")
            somthingUpdated=true
          else
            fullDetermined=false
          end 
          i=i+1
        end
      end
    end
    dump(DNIdel)
  end
end