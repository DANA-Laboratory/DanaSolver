isconstantfactor(fac::Array{Float64})=all(map(x->x==0.0,fac[2:end]))
function analysis(eqs::Equations)
  #by default expr must have a constant symbol
  push!(eqs.terms,"constant")
  j=1
  i=1
  facts::Array{Array{Float64,1},1}=Array(Array{Float64,1},0)
  for ex in eqs.exs
    fac,islinear=analysis!(ex.ex,eqs.terms,ex.termall,ex.termnonli)
    #linear equation
    if islinear 
      push!(facts,fac)
      i+=1
    else #nonlinear equations
      push!(eqs.indexnonliexs,j)
      ex.factliinnonli=copy(fac)
    end
    j+=1
  end
  vals::Array{Float64,2}=zeros(i-1,length(eqs.terms))
  for j=1:i-1 ; vals[j,1:length(facts[j])]=facts[j]; end
  eqs.terms=circshift(eqs.terms,-1)
  eqs.facts=circshift(vals,(0,-1))
  return true
end
analysis!(num::Number,terms::Array{String,1},termall::Set{String},termnonli::Set{String})=[num],true
analysis!(sym::Symbol,terms::Array{String,1},termall::Set{String},termnonli::Set{String})=analysis!(string(sym),terms,termall)
function analysis!(ssym::String,terms::Array{String,1},termall::Set{String})
  index = findfirst(terms,ssym)
  push!(termall,ssym)
  if index==0
    push!(terms,ssym)
    index=length(terms)
  end
  fac=zeros(index)
  fac[index]=1.0
  return fac,true
end
function analysis!(exp::Expr,terms::Array{String,1},termall::Set{String},termnonli::Set{String})
  if exp.head==:.
    return analysis!(string(exp),terms,termall),true
  end
  if exp.head==:call
    ln=length(exp.args)
    facAr::Array{Array{Float64,1},1}=Array(Array{Float64,1},ln-1)
    allislinear::Bool=true
    for i = [2:ln]
      facAr[i-1],islinear=analysis!(exp.args[i],terms,termall,termnonli)
      allislinear&=islinear
    end
    # fill vector with zero and same length as syms
    ln=length(terms)
    #if all(map(x->length(x)!=0,facAr))
      facAr=[vcat(fac,zeros(ln-length(fac))) for fac in facAr[1:end]]
      if exp.args[1] == :+
        sum=zeros(ln)
        for fac in facAr
          sum+=fac
        end
        return sum,allislinear
      elseif exp.args[1] == :-
        length(facAr)==1 && return (-1*facAr[1]),true
        return (facAr[1]-facAr[2]),allislinear
      elseif exp.args[1] == :*
        mult=facAr[1]
        for fac in facAr[2:end]
          if isconstantfactor(mult) == 1
            mult=mult[1]*fac
          elseif isconstantfactor(fac) == 1
            mult=mult*fac[1]
          else
            allsyms!(exp,termnonli)
            return [],false
          end
        end
        return mult,allislinear
      elseif exp.args[1] == :/
        #println(exp.args)
        #println(facAr)
        if  isconstantfactor(facAr[2]) && facAr[2][1]!=0.0
          return (facAr[1]/(facAr[2][1])),allislinear
        else
          allsyms!(exp,termnonli)
          return [],false
        end
      end
      allsyms!(exp,termnonli)
    #else
      #expr have nonlinear terms
    #end
  end
  return [],false
end
function allsyms!(exp::Expr,termnonli::Set{String})
  if exp.head == :.
    push!(termnonli,string(arg))
  elseif exp.head == :call
    for arg in exp.args[2:end]
      if typeof(arg)==Symbol
        push!(termnonli,string(arg))
      elseif typeof(arg)==Expr
        allsyms!(arg,termnonli)
      end
    end
  elseif exp.head == :tuple
    length(exp.args) == 0 && return 
    throw(DomainError())
  else
    throw(DomainError())
  end
end
