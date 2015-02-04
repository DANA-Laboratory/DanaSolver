function callfzero(fun::Function,gus::Float64,br::Vector{Float64},maxattempts=30,ord=0)
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
  try
    firstLoopDone && (return Roots.fzero(fun,gus,br))
  catch er 
    println("fzero attempt number $attempt no answer in bracket $br")
    rethrow(er)
  end
end
function callmin(eqGroup::Vector{Function},indexGroup::Vector{Tuple},lo::Vector{Float64},up::Vector{Float64},de::Vector{Float64})
  numberOfEquations=length(nonliFuns)
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
      mapreduce(x->(eqGroup[x](getindex(y,indxGroup[x])...))^2,+,[1:i])
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
