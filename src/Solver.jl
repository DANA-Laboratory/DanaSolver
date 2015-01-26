module Solver 
  export replace
  export analysis
  export replace!
  export solve
  export update!
  using DanaTypes
  using Calculus
  import DanaTypes.setfield
  #using Roots
  #replace known arguments in _eq::Function argument list and returns another function with only unknowns 
	#for unknown arguments pass NaN in _argsArray
	#getEq([1.1,2.3,NaN,-3.1,NaN],fun) -> _fun(arg1,arg2)=fun(1.1,2.3,arg1,-3.1,arg2)
	function getEq(_argsArray::Vector,_eq::Function)
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
	end

  getfield(danamodel::DanaModel,sy::Symbol)=get(Base.getfield(danamodel,sy))
  getfield(danamodel::DanaModel,ex::Expr)=getfield(getfield(danamodel,ex.args[1]),ex.args[2].value)

  setfield(danamodel::DanaModel,ex::Expr,value)=setfield(getfield(danamodel,ex.args[1]),ex.args[2].value,value)
  setfield(danamodel::DanaModel,st::String,value)=setfield(danamodel,symbol(st),value)
  
  replace(danamodel::DanaModel,num::Float64)=num
  replace(danamodel::DanaModel,num::Int32)=num
  replace(danamodel::DanaModel,num::Int64)=num

  function replace(danamodel::DanaModel)
    new_equations::Array{Expr,1}=Array(Expr,0)
    for exp in danamodel.equationsFlow
      push!(new_equations,replace(danamodel,exp))
    end
    return new_equations
  end  
  
  function replace(danamodel::DanaModel,sym::Symbol)
    var = getfield(danamodel,sym)
    return (isunknown(var) ? sym : var)
  end
  #replace symbols with their values & := with :-
  function replace(danamodel::DanaModel,exp::Expr)
    if exp.head == :call
      new_exp=Expr(:call,exp.args[1])
      for expr in exp.args[2:end]
        push!(new_exp.args,replace(danamodel,expr))
      end
      return new_exp
    elseif exp.head == :(=)
      return Expr(:call,:-,replace(danamodel,exp.args[1]),replace(danamodel,exp.args[2]))
    elseif exp.head == :.
      var = getfield(danamodel,exp)
      return (isunknown(var) ? exp : var)
    end
    return exp
  end

  function analysis(exps::Array{Expr,1})
    syms::Array{String,1}=Array(String,0)
    nolinearFunctions::Array{Function,1}=Array(Function,0)
    nolinearArgs::Array{Set{String},1}=Array(Set{String},0)
    #by default expr must have a constant symbol
    push!(syms,"constant")
    i=1
    facs::Array{Array{Float64,1},1}=Array(Array{Float64,1},length(exps))
    for exp in exps
      symsLoc::Set{String}=Set{String}()
      fac=analysis!(exp,syms,symsLoc)
      #linear equation
      if length(fac)>0
        facs[i]=fac
        i+=1
      else #nonlinear
        nonlinearExpr::Expr=:(()->())
        nonlinearExpr.args[2].args[2]=exp
        for s in symsLoc
          push!(nonlinearExpr.args[1].args,symbol(s))
        end
        push!(nolinearFunctions, eval(nonlinearExpr))
				push!(nolinearArgs, symsLoc)
      end
    end
    vals::Array{Float64,2}=zeros(i-1,length(syms))
    for j=1:i-1 ; vals[j,1:length(facs[j])]=facs[j]; end
    return  circshift(vals,(0,-1)) , circshift(syms,-1) , nolinearFunctions , nolinearArgs
  end
  
  analysis!(num::Number,syms::Array{String,1},symsLoc::Set{String})=[num]
  
  analysis!(sym::Symbol,syms::Array{String,1},symsLoc::Set{String})=analysis!(string(sym),syms,symsLoc)
  
  function analysis!(ssym::String,syms::Array{String,1},symsLoc::Set{String})
    index = findfirst(syms,ssym)
    push!(symsLoc,ssym)
    if index==0
      push!(syms,ssym)
      index=length(syms)
    end
    fac=zeros(index)
    fac[index]=1.0
    return fac
  end
  
  function analysis!(exp::Expr,syms::Array{String,1},symsLoc::Set{String})
    if exp.head==:.
      return analysis!(string(exp),syms,symsLoc)
    end
    if exp.head==:call
      ln=length(exp.args)
      facAr::Array{Array{Float64,1},1}=Array(Array{Float64,1},ln-1)
      for i = [2:ln]
        facAr[i-1]=analysis!(exp.args[i],syms,symsLoc)
      end
      # fill vector with zero and same length as syms
      ln=length(syms)
      if all(map(x->length(x)!=0,facAr))
        facAr=[vcat(fac,zeros(ln-length(fac))) for fac in facAr[1:end]]
        if exp.args[1] == :+ 
          sum=zeros(ln)
          for fac in facAr
            sum+=fac
          end
          return sum
        elseif exp.args[1] == :-
          return (facAr[1]-facAr[2])
        elseif exp.args[1] == :*
          mult=facAr[1]
          for fac in facAr[2:end]
            if isconstantfactor(mult) == 1 
              mult=mult[1]*fac
            elseif isconstantfactor(fac) == 1 
              mult=mult*fac[1]
            else
              return []
            end
          end
          return mult
        elseif exp.args[1] == :/
          if  isconstantfactor(facAr[2])
            return (facAr[1]/(facAr[2][1]))
          else
            return []
          end
        end
      end
    end
    return []
  end
  
  isconstantfactor(fac::Array{Float64})=all(map(x->x==0.0,fac[2:end]))

  #solve equations of a model
  #call solve after setEquationFlow()
  function solve(danamodel::DanaModel)
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
        end
      catch y
        if isa(y,DomainError)
          println ("can't simplify following equation:\n",eq);
        end
        println ("can't simplify following equation:\n",eq);
        throw(y)
      end
    end
    #equations=vals*vars (linear equations)
    vals,vars,nolinearFunctions,nolinearArgs=analysis(sequations)
    #reduced row echelon form
    rreModel=rref(vals)
    return rreModel,vars,nolinearFunctions,nolinearArgs
  end

  function update!(danamodel::DanaModel,rre::Array{Float64,2},vars::Array{String,1})
    somthingUpdated=false
		fullDetermined=true
    si=size(rre)
    clms=si[2]
    rows=si[1]
		if clms<2 || rows<1
			return somthingUpdated,false
		end
		if rows!=clms-1
			fullDetermined=false
		end
    valz=slicedim(rre,2,clms)
    lda=1
    for r = 1:rows
      arow=slice(rre,r,1:clms-1)
      while lda<clms && arow[lda]==0
        lda+=1
      end
      j=lda+1
      while j<clms && arow[j]==0
        j+=1
      end
      if j==clms
        setfield(danamodel,vars[lda],-1*valz[r])
        somthingUpdated=true
      else
        fullDetermined=false
      end
    end
    return somthingUpdated,fullDetermined
  end
  
  #replace variables in an expresion with floats
  function replace!(expr::Expr,m::Dict{String,Float64})
    strerp::String=string(expr)
    if haskey(m,strerp)
      expr=m[strerp]
    elseif isa(expr,Expr)
      for i = 1:length(expr.args)
        strerp=string(expr.args[i])
        if haskey(m,strerp)
          expr.args[i] =  m[strerp]
        elseif isa(expr.args[i],Expr) && expr.args[i].head!=:.
          replace!(expr.args[i],m)
        end
      end
    end
  end

  #replace variables in expresion with floats
  function replace!(expr::Expr,exprVars::Dict{Expr,Float64},symVars::Dict{Symbol,Float64})
    if haskey(exprVars,expr)
      expr=exprVars[strerp]
    else isa(expr,Expr)
      for i = 1:length(expr.args)
        if isa(expr.args[i],Expr) && haskey(exprVars,expr.args[i])
          expr.args[i] =  exprVars[expr.args[i]]
        elseif isa(expr.args[i],Expr) && expr.args[i].head!=:.
          replace!(expr.args[i],exprVars,symVars)
        elseif isa(expr.args[i],Symbol) && haskey(symVars,expr.args[i])
          expr.args[i] =  symVars[expr.args[i]]
        end
      end
    end
  end
  
  #replace a variable symbol in expresion with float
  function replace!(expr::Expr,var::Symbol,val::Float64)
    for i = 1:length(expr.args)
      if isa(expr.args[i],Symbol) && expr.args[i]==var
        expr.args[i] =  val
      elseif isa(expr.args[i],Expr) && expr.args[i].head!=:.
        replace!(expr.args[i],var,val)
      end
    end
  end

  #replace a variable expression in a expression with float
  function replace!(expr::Expr,var::Expr,val::Float64)
    if expr==var
      expr=var
    else
      for i = 1:length(expr.args)
        if isa(expr.args[i],Expr)
          if expr.args[i]==var
            expr.args[i] =  val
          elseif expr.args[i].head!=:.
            replace!(expr.args[i],var,val)
          end
        end
      end
    end
  end
end
