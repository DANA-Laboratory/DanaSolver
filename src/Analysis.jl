isconstantfactor(fac::Array{Float64})=all(map(x->x==0.0,fac[2:end]))
function analysis(exps::Array{Expr,1})
  syms::Array{String,1}=Array(String,0)
  nonlinearExpIndx::Array{Int,1}=Array(Int,0)
  symsLocArray::Array{Set{String},1}=Array(Set{String},0)
  symsNonlinearArray::Array{Set{String},1}=Array(Set{String},0)
  #by default expr must have a constant symbol
  push!(syms,"constant")
  i=1
  j=1
  facs::Array{Array{Float64,1},1}=Array(Array{Float64,1},length(exps))
  nonlifacs::Array{Array{Float64,1},1}=Array(Array{Float64,1},length(exps))
  for exp in exps
    symsLoc::Set{String}=Set{String}()
    symsNonlinear::Set{String}=Set{String}()
    fac,islinear=analysis!(exp,syms,symsLoc,symsNonlinear)
    push!(symsLocArray, symsLoc)
    push!(symsNonlinearArray, symsNonlinear)
    #linear equation
    if islinear 
      facs[i]=fac
      i+=1
    else #nonlinear equations
      push!(nonlinearExpIndx,j)
      nonlifacs[j]=fac
    end
    j+=1
  end
  vals::Array{Float64,2}=zeros(i-1,length(syms))
  for j=1:i-1 ; vals[j,1:length(facs[j])]=facs[j]; end
  return  circshift(vals,(0,-1)) , circshift(syms,-1) , nonlinearExpIndx , symsLocArray , symsNonlinearArray , nonlifacs
end
analysis!(num::Number,syms::Array{String,1},symsLoc::Set{String},symsNonlinear::Set{String})=[num],true
analysis!(sym::Symbol,syms::Array{String,1},symsLoc::Set{String},symsNonlinear::Set{String})=analysis!(string(sym),syms,symsLoc)
function analysis!(ssym::String,syms::Array{String,1},symsLoc::Set{String})
  index = findfirst(syms,ssym)
  push!(symsLoc,ssym)
  if index==0
    push!(syms,ssym)
    index=length(syms)
  end
  fac=zeros(index)
  fac[index]=1.0
  return fac,true
end
function analysis!(exp::Expr,syms::Array{String,1},symsLoc::Set{String},symsNonlinear::Set{String})
  if exp.head==:.
    return analysis!(string(exp),syms,symsLoc),true
  end
  if exp.head==:call
    ln=length(exp.args)
    facAr::Array{Array{Float64,1},1}=Array(Array{Float64,1},ln-1)
    allislinear::Bool=true
    for i = [2:ln]
      facAr[i-1],islinear=analysis!(exp.args[i],syms,symsLoc,symsNonlinear)
      allislinear&=islinear
    end
    # fill vector with zero and same length as syms
    ln=length(syms)
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
            allsyms!(exp,symsNonlinear)
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
          allsyms!(exp,symsNonlinear)
          return [],false
        end
      end
      allsyms!(exp,symsNonlinear)
    #else
      #expr have nonlinear terms
    #end
  end
  return [],false
end
function allsyms!(exp::Expr,symsNonlinear::Set{String})
  if exp.head == :.
    push!(symsNonlinear,string(arg))
  elseif exp.head == :call
    for arg in exp.args[2:end]
      if typeof(arg)==Symbol
        push!(symsNonlinear,string(arg))
      elseif typeof(arg)==Expr
        allsyms!(arg,symsNonlinear)
      end
    end
  elseif exp.head == :tuple
    length(exp.args) == 0 && return 
    throw(DomainError())
  else
    throw(DomainError())
  end
end
