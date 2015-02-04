
setfield!(danamodel::DanaModel,ex::Expr,value)=setfield!(getfield(danamodel,ex.args[1]),ex.args[2].value,value)
setfield!(danamodel::DanaModel,st::String,value)=setfield!(danamodel,symbol(st),value)

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
      try
        setfield!(danamodel,vars[lda],-1*valz[r])
        somthingUpdated=true
      catch er
        println("in update! can't set ",vars[lda]," with ",-1*valz[r])
        throw(er)
      end
    else
      fullDetermined=false
    end
  end
  return somthingUpdated,fullDetermined
end
