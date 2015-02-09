isconstantfactor(fac::Array{Float64})=all(map(x->x==0.0,fac[2:end]))
function exprTofunction(ex::Expr,symsLoc::Set{String})
  ret::Expr=:(()->())
  ret.args[2].args[2]=ex
  for s in symsLoc
    push!(ret.args[1].args,symbol(s))
   end
  return eval(ret)
end
function analysis(exps::Array{Expr,1},all::Bool=false)
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
    #nonlinear or all
    if all || length(fac)==0
      push!(nolinearFunctions, exprTofunction(exp,symsLoc))
      push!(nolinearArgs, symsLoc)
    end
    #linear equation
    if length(fac) != 0
      facs[i]=fac
      i+=1
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
        length(facAr)==1 && return (-1*facAr[1])
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
