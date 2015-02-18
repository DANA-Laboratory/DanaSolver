const MaxAttempt=100
initial=NaN 
function callffzero(f::Function,finitial::Float64,fbracket::Vector{Float64},g::Function,ginitial::Float64,gbracket::Vector{Float64})
  global initial
  initial=NaN
  re=callfzero(x->getfxy0whengzeroatx(f,g,x,ginitial,gbracket),finitial,fbracket,500)
  return re,initial
end
#find y0 when y->g(x,y)==0 and return f(x,y0)
function getfxy0whengzeroatx(f::Function,g::Function,x::Float64,ginitial,gbracket)
  attempt=0
  global initial
  isnan(initial) && (initial=ginitial)
  while attempt<MaxAttempt
    try
      y0=Roots.fzero(y->g([x,y]),initial,order=0)
      #println("find zero @(x=$x,y=$y0) initial=$initial f([x,y])=",f([x,y0]))
      initial=y0
      return f([x,y0])
    catch er
      if er == DomainError()
        initial=rand()*(gbracket[2]-gbracket[1])+gbracket[1]
        attempt+=1
      else
        rethrow(er)
      end
    end
  end
  #println ("in callffzero no answer after $attempt attempts")
  throw(DomainError)
end
function callfzero(fun::Function,gus::Float64,br::Vector{Float64},maxattempts=MaxAttempt,ord=0)
  result=NaN
  attempt=0
  firstLoopDone=false
  while (!firstLoopDone && attempt<maxattempts)
    attempt+=1
    try
      result=Roots.fzero(fun,gus,order=ord)
      if  (!isnan(result) && result<=br[2] && result>=br[1])
        return result,attempt
      else
        firstLoopDone=true # no more chanse
      end
    catch er
      if  attempt==maxattempts
        println("fzero attempt number $attempt fail result=$result gus=$gus")
        rethrow(er)
      else
        gus=rand()*(br[2]-br[1])+br[1]
      end
    end
  end
  #println("in callfzero attemp=$attempt firstLoopDone=$firstLoopDone result=$result")
  try
    firstLoopDone && (return Roots.fzero(fun,gus,br))
  catch er 
    println("fzero attempt number $attempt no answer in bracket $br")
    rethrow(er)
  end
end
function callmin(eqGroup::Vector{Function},indxGroup::Vector{Vector{Int}},lo::Vector{Float64},up::Vector{Float64},de::Vector{Float64})
  numberOfEquations=length(eqGroup)
  somethingUpdated=false
  opt = Opt(:GN_DIRECT_L, length(eqGroup))
  lower_bounds!(opt, lo)
  upper_bounds!(opt, up)
  stopval!(opt, 1.0e-12)
  maxtime!(opt, 1.0*length(eqGroup))
  ftol_abs!(opt, 1.0e-19)
  ftol_rel!(opt, 1.0e-18)
  optfun=(y,gradient)->begin
    try
      mapreduce(x->(eqGroup[x](getindex(y,indxGroup[x])...))^2,+,[1:numberOfEquations])
    catch er
      println("in nonlinear optimize , fail with following vals: ",y)
      rethrow(er)
    end
  end
  min_objective!(opt,optfun)
  (minf,minx,ret)=optimize(opt,de)
  if "$ret"=="STOPVAL_REACHED"
    return minx
  elseif "$ret"=="MAXTIME_REACHED"
    println("NLopt fail to MAXTIME_REACHED")
    return nothing
  end
end
