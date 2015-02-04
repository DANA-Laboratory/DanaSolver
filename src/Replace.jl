
replace(danamodel::DanaModel,num::Float64)=num
replace(danamodel::DanaModel,num::Int32)=num
replace(danamodel::DanaModel,num::Int64)=num

getfield(danamodel::DanaModel,sy::Symbol)=get(Base.getfield(danamodel,sy))
getfield(danamodel::DanaModel,ex::Expr)=getfield(getfield(danamodel,ex.args[1]),ex.args[2].value)
    
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

# replace all occurences of var with wexpr 
function replace!(ex:Expr,var::Symbol,wexpr::Expr)
  le=length(ex.args)
  i=1
  while(i<le)
    if (typeof(ex.args[i])==Symbol && ex.args[i]==var)
      ex.args[i]=wexpr
    elseif (typeof(ex.args[i])==Expr)
      replace!(ex.args[i],var,wexpr)
    end
    i+=1
  end
end
