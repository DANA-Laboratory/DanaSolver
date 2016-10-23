module Solver
  using DanaTypes
  using Calculus
  using Roots
  using NLopt
  export solve!
  import DanaTypes.setfield!
  include("Replace.jl")
  include("Equations.jl")
  include("Analysis.jl")
  include("Update.jl")
  include("Calls.jl")

  #convert expr to function of symsLoc --moves from analysis
  function exprTofunction(ex::Expr,symsLoc::Vector{String})
    ret::Expr=:(()->())
    ret.args[2].args[2]=ex
    for s in symsLoc
      push!(ret.args[1].args,symbol(s))
    end
    return eval(ret)
  end

  #replace known arguments in _eq::Function argument list and returns another function with only unknowns
  #for unknown arguments pass NaN in _argsArray
  #getEq([1.1,2.3,NaN,-3.1,NaN],fun) -> _fun(arg1,arg2)=fun(1.1,2.3,arg1,-3.1,arg2)

  #=	function getEq(_argsArray::Vector,_eq::Function)
    ex=:((argsArray)->$_eq())
    exx=[_eq]
    len=length(_argsArray)
    j=1
    for i in 1:len
      if isnan(_argsArray[i])
        exx=vcat(exx,:(argsArray[$j]))
        j+=1
      else
        exx=vcat(exx,_argsArray[i])
      end
    end
    ex.args[2].args[2].args=exx
    return eval(ex)
  end =#

  #solve li plus one nonlinear
  function sliponl!(danamodel::DanaModel)
    somethingUpdated=true
    fullDetermined=false
    noliTrys=0
    nonlTrys=0
    df_gus=1.0
    df_br=[0.0,2.0e10]
    while (true)
      somethingUpdated,fullDetermined,eqs,noTrys=slstsubnfd!(danamodel)
      noliTrys+=noTrys
      if(!fullDetermined)
        #find a system to solve
        args=Array(Set{String},0)
        #for ex in eqs.exs; push!(args,ex.termall); end;
        varIndex,allVars,eqIndex=findsystem(eqs,1)
        if (varIndex!=false)
          unknown=allVars[varIndex[1]]
          println(unknown)
          println(eqIndex)
          println(allVars)
          println(varIndex)
          fun=Solver.exprTofunction(eqs.exs[eqIndex[1]].ex,unknown)
          result,attempt=callfzero(fun,getdefault(danamodel,unknown[1]),getbracket(danamodel,unknown[1]))
          Solver.setfield!(danamodel,unknown[1],result)
          setEquationFlow(danamodel)
          somethingUpdated=true
          nonlTrys+=1
        else
          varIndex,allVars,eqIndex=findsystem(eqs,2)
          if (varIndex!=false)
            #fails to find a nonlinear equation of one unknown
            #this system can simplified to a nonlinear equation of one unknown
            #println(equations)
            println("eqIndex=",eqIndex)
            for ind in eqIndex
              println(eqs.exs[ind].ex)
              println(eqs.exs[ind].termnonli," ",eqs.exs[ind].factliinnonli)
            end
            println("here")
          end
          somethingUpdated=false
        end
      end
      if (!somethingUpdated || fullDetermined)
        return somethingUpdated,fullDetermined,noliTrys,nonlTrys
      end
    end
  end

  #main loop
  function solve!(danamodel::DanaModel)
    somethingUpdated=true
    fullDetermined=false
    nonlTrys=0
    noliTrys=0
    while (somethingUpdated && !fullDetermined)
      somethingUpdated,fullDetermined,noliTrys,nonlTrys=sliponl!(danamodel)
      if !fullDetermined
        #somethingUpdated,fullDetermined=ssonle!(danamodel,nonliFuns,nonliArgs)
      end
    end
    return somethingUpdated,fullDetermined,noliTrys,nonlTrys
  end

  # solve system of nonlinear equations
  function ssonle!(danamodel::DanaModel,nonliFuns::Array{Function,1},nonliArgs::Array{Set{String},1})
    somthingUpdated=false
    fullDetermined=false
    numberOfEquations=length(nonliFuns)
    noe=2 #number of equations
    while (noe<=numberOfEquations && !somthingUpdated)
      eqIndexes=getapsoe(1,numberOfEquations,noe)
      for eqIndex in eqIndexes
        varGroup=getindex(nonliArgs,eqIndex)
        allVars=union(varGroup...)
        if length(allVars) == noe
          eqGroup=getindex(nonliFuns,eqIndex)
          #println("eqGroup=",eqIndex," for vars:",allVars)
          indxGroup=map(x->indexin([x...],[allVars...]),varGroup)
          opt = Opt(:GN_DIRECT_L, noe)
          lower_bounds!(opt, [1.0e-3, 1.0])
          upper_bounds!(opt, [10.0,2500])
          stopval!(opt, 1.0e-12)
          maxtime!(opt, 1.0*noe)
          #ftol_abs!(opt, 1.0e-19)
          #ftol_rel!(opt, 1.0e-18)
          min_objective!(opt, (y,gradient)->mapreduce(x->(apply(eqGroup[x],getindex(y,indxGroup[x])))^2,+,[1:noe]))
          (minf,minx,ret)=optimize(opt,ones(Float64,noe))
          #println("got $minf at $minx (returned $ret)")
          if "$ret"=="STOPVAL_REACHED"
            for j in [1:noe]
              setfield!(PR,[nonliVars[eqIndex[1]]...][j],minx[j])
            end
            somthingUpdated=true
            if noe==numberOfEquations
              fullDetermined=true
            end
            break
          end
        end
      end
    end
    return somthingUpdated,fullDetermined
  end

  # Solve Linear System Til Something Updated But Not Full Determined
  # setEquationFlow -> loop(replace->symplify->analysis->update)
  function slstsubnfd!(danamodel::DanaModel)
    noTrys=0
    somethingUpdated=true
    fullDetermined=false
    setEquationFlow(danamodel)
    while (true)
      noTrys+=1
      rVls,eqs=solvelinear(danamodel)
      somethingUpdated,fullDetermined=update!(danamodel,rVls,eqs.terms)
      somethingUpdated && setEquationFlow(danamodel)
      if (!somethingUpdated || fullDetermined)
        return somethingUpdated,fullDetermined,eqs,noTrys
      end
    end
    #return somethingUpdated,fullDetermined,nonliExprIndx,args,equations,noTrys,symsNonlinearArray,nonlifacs
  end

  #solve equations of a model
  #call solve after setEquationFlow()
  function solvelinear(danamodel::DanaModel)
    #replace variables & parameters with values
    equations=replace(danamodel)
    #call simplify on eatch equation
    sequations::Array{Expr,1}=Array(Expr,0)
    for eq in equations
      try
        seq=simplify(eq)
        if isa(seq,Expr)
          push!(sequations,seq)
        else
          # nothing to do!
        end
      catch y
        if isa(y,DomainError)
          println("can't simplify following equation:\n",eq);
        end
        println("can't simplify following equation:\n",eq);
        throw(y)
      end
    end
    #equations=vals*vars (linear equations)
    eqs=Equations(map(Equation,sequations))
    analysis(eqs)
    #reduced row echelon form
    rreModel=rref(eqs.facts)
    return rreModel,eqs
  end
  function rref(U::Array{Float64,2})
    nr, nc = size(U)
    e = eps(norm(U,Inf))
    i = j = 1
    while i <= nr && j <= nc
      (m, mi) = findmax(abs(U[i:nr,j]))
      mi = mi+i - 1
      if m <= e
        U[i:nr,j] = 0
        j += 1
      else
        for k=j:nc
          U[i, k], U[mi, k] = U[mi, k], U[i, k]
        end
        d = U[i,j]
        for k = j:nc
          U[i,k] /= d
        end
        for k = 1:nr
          if k != i
            d = U[k,j]
            for l = j:nc
              U[k,l] -= d*U[i,l]
            end
          end
        end
        i += 1
        j += 1
      end
    end
    U
  end
  rref(x::Number) = one(x)

  function findsystem(eqs::Equations,nonliTerms=1)
    numberOfEquations=length(eqs.exs)
    noe=nonliTerms #number of equations
    while (noe<=numberOfEquations)
      eqIndexes=getapsoe(1,numberOfEquations,noe)
      for eqIndex in eqIndexes
        varGroup=Array(Set{String},0)
        for ind in eqIndex ; push!(varGroup,eqs.exs[ind].termnonli) ; end
        allVars=union(varGroup...)
        if (length(allVars) == nonliTerms)
          varGroup=Array(Set{String},0)
          for ind in eqIndex ; push!(varGroup,eqs.exs[ind].termall) ; end
          allVars=[union(varGroup...)...]
          if (length(allVars) == noe)
            if(noe>nonliTerms)#remove=noe-nonliTerms
              way=Dict{String,Array{Int,1}}{}
            end
            varIndex=map(x->indexin([x...],allVars),varGroup)
            return varIndex,allVars,eqIndex
          end
        end
      end
      noe+=1
    end
    return false,false,false
  end

  #generate indexes for all possible system of _noe equations[APSOE].
  #where equations are selected from _minIndex to _maxIndex of a list of equations
  #TODO use combinators
  function getapsoe(minIndex::Int,maxIndex::Int,noe::Int)
    if 1<noe
      jj::Vector=Vector[]
      for k in [minIndex+1:maxIndex]
        j=getapsoe(k,maxIndex,noe-1)
        jj=append!(jj,[push!(e,k-1) for e in j])
      end
      return jj
    else
      return [[i] for i in minIndex:maxIndex]
    end
  end

  #solve an expr for var with factor->fac
  #must be linear in term of var
  function solveexpr!(expr::Expr,var::Symbol,fac::Float64)
    replace!(expr,var,0.0)
    expr=simplify(expr)
    expr=Expr(:call,:/,:(-1*$expr),fac)
  end

  function whatifcombineeqs(lia::Set{String},nlia::Set{String},lib::Set{String},nlib::Set{String})
    nli_re=union(nlia,nlib)
    com=[intersect(lia,lib),intersect(lia,nlib),intersect(nlia,lib)]
    un=[[union(lia,lib) union(nlia,nlib)];[lib union(lia,nlib,nlia)];[lia union(nlia,lib,nlib)]]
    #number of possible scenarions=length(common(lia,lib))
    re::Array{Set{String},1}=Array(Set{String},0)
    j=1
    for i in [1:3]
      for c in (com[i])
        lic=copy(un[i,1])
        nlic=copy(un[i,2])
        delete!(lic,c)
        delete!(nlic,c)
        setdiff!(lic,nlic)
        push!(re,lic,nlic)
      end
    end
    return re
  end
end
