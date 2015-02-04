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
