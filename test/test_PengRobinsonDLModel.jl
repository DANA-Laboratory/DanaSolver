# REF[1] Chemical Process Design and Integration By Robin Smith
# REF[2] Engineering and Chemical Thermodynamics By Milo D. Koretsky
# REF[3] Perry HandBook

#generate indexes for all possible system of _noe equations[APSOE].
#where equations is a selection from _minIndex to _maxIndex of a list of equations
using NLopt
function dumpme(var)
  println("****** dump ******")
  println("PR.Tc , PR.Pc , PR.af , PR.Zc , PR.h_Dep , PR.s_Dep , PR.g_Dep , PR.Pr , PR.Tr , PR.vr")
  println(get(var.Tc)," , ",get(var.Pc)," , ",get(var.af)," , ",get(var.Zc)," , ",get(var.h_Dep)," ,  ",get(var.s_Dep)," , ",get(var.g_Dep)," , ",get(var.Pr)," , ",get(var.Tr)," , ",get(var.vr))
end
function testPRDLvriousunknowns()
  println("**** test RP DL variouns unknowns ****")
  nonliFuns::Array{Function,1}=Array(Function,0)
  nonliVars::Array{Set{String},1}=Array(Set{String},0)
  for k in [1:3]
    cNo="75-07-0" #Acetaldehyde
    PR=DANAPengRobinsonDL()
    Tc,Pc,af=getvalueforname("Criticals","Acetaldehyde")
    setfield!(PR,:Tc,Tc)
    setfield!(PR,:Pc,Pc)
    setfield!(PR,:af,af)
    h_Dep,s_Dep,Pr,Tr,vr=[(NaN,NaN,1.0,1.0,NaN),(NaN,NaN,NaN,1.0,1.0),(NaN,NaN,1.0,NaN,1.0),(1.1096883953196783e7,NaN,1.0,NaN,NaN),(1.1096883953196783e7,NaN,NaN,1.0,NaN),(1.1096883953196783e7,NaN,NaN,NaN,0.21722233067387567),(NaN,20135.0,1.0,NaN,NaN),(NaN,20135.0,NaN,1.0,NaN),(NaN,20135.0,NaN,NaN,0.21722233067387567)][k]
    setfield!(PR,:h_Dep,h_Dep)
    setfield!(PR,:s_Dep,s_Dep)
    setfield!(PR,:Pr,Pr)
    setfield!(PR,:Tr,Tr)
    setfield!(PR,:vr,vr)
    println("----------Solution Starts----------")
    dumpme(PR)
    somethingUpdated,fullDetermined,noliTrys,nonlTrys = solve!(PR)
    println("--somethingUpdated=$somethingUpdated,fullDetermined=$fullDetermined,noliTrys=$noliTrys,nonlTrys=$nonlTrys--")
    dumpme(PR)
  end
end
#=
    somethingUpdated=true
    fullDetermined=false
    while (somethingUpdated && !fullDetermined)
      while (somethingUpdated && !fullDetermined)
        setEquationFlow(PR);
        #linear solver
        somethingUpdated,fullDetermined,nonliFuns,nonliVars,noTrys=Solver.slstsubnfd!(PR)
        println("Linear Solver. number of attempts=",noTrys)
        if !fullDetermined
          somethingUpdated=false
          i=1
          fullDetermined=true
          #search for non-linear equations with only one unknown
          numberOfEquations=length(nonliFuns)
          while (i<=numberOfEquations && !somethingUpdated)
            if length(nonliVars[i])==1
              f_res=NaN
              varName=[nonliVars[i]...][1]
              br=getbracket(PR,varName)
              attem=0
              gus=getdefault(PR,varName)
              result=NaN
              while (isnan(result) || result>br[2] || result<br[1])
                try
                  attem+=1
                  result=Roots.fzero(nonliFuns[i],gus,order=0)
                  f_res=nonliFuns[i](result)
                  if (!isnan(result) && result<=br[2] || result>=br[1])
                    somethingUpdated=true
                    println("fzero attempt number $attem succ result=$result gus=$gus")
                  else
                    println("fzero attempt number $attem fail result=$result gus=$gus")
                    gus=rand()*(br[2]-br[1])+br[1]
                  end
                catch
                  gus=rand()*(br[2]-br[1])+br[1]
                end
              end
              if somethingUpdated==true
                Solver.setfield!(PR,varName,result)
                setEquationFlow(PR);
                println("non-Linear Solver, solves one equation. $varName = $result  f(result)=$f_res")
              else
                println("fzero fail")
              end
            else
              fullDetermined=false
            end
            i=i+1
          end
        end
      end
      if !fullDetermined
        numberOfEquations=length(nonliFuns)
        println("non-Linear multi-Equation Solver.")
        #fail to fined non-linear equations with only one unknown
        somethingUpdated=false
        i=2
        println(nonliVars)
        while (i<=numberOfEquations && !somethingUpdated)
          println("try system of ",i," equations")
          eqIndexes=getAPSOE(1,numberOfEquations,i)
          for eqIndex in eqIndexes
            varGroup=getindex(nonliVars,eqIndex)
            allVars=union(varGroup...)
            i==2 && (println(varGroup,' ',length(allVars),' ',eqIndex))
            if length(allVars) == i
              eqGroup=getindex(nonliFuns,eqIndex)
              println("eqGroup=",eqIndex," for vars:",allVars)
              lo=Array(Float64,0)
              up=Array(Float64,0)
              de=Array(Float64,0)
              for ii in 1:i
                push!(lo,(getbracket(PR,[allVars...][ii]))[1])
                push!(up,(getbracket(PR,[allVars...][ii]))[2])
                push!(de,(getdefault(PR,[allVars...][ii])))
              end
              indxGroup=map(x->indexin([x...],[allVars...]),varGroup)
              println("lower bounds=",lo)
              println("upper bounds=",up)
              println("defaults=",de)
              println("indexGroup=",indxGroup)
              # use multiple fzeroz instead
              if (i==2)
                funcs=[y->eqGroup[1](getindex(y,indxGroup[1])...),y->eqGroup[2](getindex(y,indxGroup[2])...)]
                println("using multiple fzeroz")
                @time r=Solver.callfzero(x->Solver.callffzero(x,funcs[1],de[2],[lo[2],up[2]],funcs[2]),de[1],[lo[1],up[1]],1000)
                println("result=",r[1])
                println("using opt")
                opt = Opt(:GN_DIRECT_L, i)
                lower_bounds!(opt, lo)
                upper_bounds!(opt, up)
                stopval!(opt, 1.0e-12)
                maxtime!(opt, 1.0*i)
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
                @time(minf,minx,ret)=optimize(opt,de)
                println("got $minf at $minx(returned $ret)")
                if "$ret"=="STOPVAL_REACHED"
                  for j in [1:i]
                    Solver.setfield!(PR,[nonliVars[eqIndex[1]]...][j],minx[j])
                  end
                  somethingUpdated=true
                  if i==numberOfEquations
                    fullDetermined=true
                  end
                  break
                elseif "$ret"=="MAXTIME_REACHED"
                  println("NLopt fail to MAXTIME_REACHED")
                  return nothing
                end
              end
            end
          end
          i+=1
        end
      end
    end
    if fullDetermined
      dumpme(PR)
      println("----------Solution Done----------")
    end
  end
end
#*********************
function testVariousKnowns()
  #P & T
  # butane # 106-97-8
  cNo="75-07-0"
  v=0.0;
  PR=DANAPengRobinson()
  PR.Tc,PR.Pc,PR.af=getValueForCasNo("Criticals",cNo)
  PR.P=PR.Pc
  PR.T=PR.Tc
  somethingUpdated=true
  fullDetermined=false
  while (somethingUpdated && !fullDetermined)
    setEquationFlow(PR);
    rVls,vars,nonliFuns,nonliVars=solve(PR)
    somethingUpdated,fullDetermined=update!(PR,rVls,vars)
    if fullDetermined
      println("solved for v! PR.v=",PR.v)
      v=PR.v
    end
  end
  #v & T
  PR=DANAPengRobinson()
  PR.Tc,PR.Pc,PR.af=getValueForCasNo("Criticals",cNo)
  PR.v=v
  PR.T=PR.Tc
  somethingUpdated=true
  fullDetermined=false
  while (somethingUpdated && !fullDetermined)
    setEquationFlow(PR);
    rVls,vars,nonliFuns,nonliVars=solve(PR)
    somethingUpdated,fullDetermined=update!(PR,rVls,vars)
    if fullDetermined
      println("solved for P! PR.P=",PR.P)
    end
  end
  #v & P
  PR=DANAPengRobinson()
  PR.Tc,PR.Pc,PR.af=getValueForCasNo("Criticals",cNo)
  PR.v=v
  PR.P=PR.Pc
  somethingUpdated=true
  fullDetermined=false
  while (somethingUpdated && !fullDetermined)
    setEquationFlow(PR);
    rVls,vars,nonliFuns,nonliVars=solve(PR)
    somethingUpdated,fullDetermined=update!(PR,rVls,vars)
    if !fullDetermined
      i=1
      fullDetermined=true
      while (i<=length(nonliFuns))
        if length(nonliVars[i])==1
          result=Roots.fzero(nonliFuns[i],[0,typemax(Int64)])
          HelperEquation.setfield!(PR,nonliVars[i][1],result)
          somethingUpdated=true
        else
          fullDetermined=false
        end
        i=i+1
      end
    end
    if fullDetermined
      println("solved for T! PR.T=",PR.T)
    end
  end
  #h & P
  global y=PR
  res=optimize(optFunctionHP, [PR.T/20,PR.h_Dep/20]);
  return res;
end
function optFunctionHP(x::Vector)
  b=y.b
  R=y.R
  Tc=y.Tc
  v=y.v
  d=y.d
  k=y.k
  P=y.P
  T=x[1];
  h_Dep=x[2];
  ret1=h_Dep-((-4*(b^3*R*T*Tc-2*b^2*R*T*Tc*v+d*(Tc-2*k*(-1+sqrt(T/Tc))*Tc+k^2*(T+Tc-2*sqrt(T/Tc)*Tc))*v^2-b*v*(d*(Tc-2*k*(-1+sqrt(T/Tc))*Tc+k^2*(T+Tc-2*sqrt(T/Tc)*Tc))+R*T*Tc*v)))/(Tc*(b-v)*(b^2-2*b*v-v^2))-(sqrt(2)*d*(1+k)*(-1+k*(-1+sqrt(T/Tc)))*log(-1+(b+v)/(sqrt(2)*b)))/b+(sqrt(2)*d*(1+k)*(-1+k*(-1+sqrt(T/Tc)))*log(1+(b+v)/(sqrt(2)*b)))/b)/4;
  ret2=P-R*T/(v-b)+(d*(1+k*(1-sqrt(T/Tc)))^2)/(v*v+2*b*v-b*b);
  return ret1^2+ret2^2;
end
# Verification REF[2] Example 5.4
function testDeparture()
  DNpr1=DANAPengRobinson()
  DNpr2=DANAPengRobinson()
  DNpr1.P=9.47*1e5
  DNpr1.T=80+273.15
  DNpr2.P=18.9*1e5
  DNpr2.T=120+273.15
  # butane
  tc,pc,af=getValueForCasNo("Criticals","106-97-8")
  DNpr1.Tc=tc
  DNpr1.Pc=pc
  DNpr1.af=af
  DNpr2.Tc=tc
  DNpr2.Pc=pc
  DNpr2.af=af
  somethingUpdated=true
  fullDetermined=false
  while (somethingUpdated && !fullDetermined)
    setEquationFlow(DNpr1);
    rVls,vars,nonliFuns,nonliVars=solve(DNpr1)
    somethingUpdated,fullDetermined=update!(DNpr1,rVls,vars)
  end
  somethingUpdated=true
  fullDetermined=false
  while (somethingUpdated && !fullDetermined)
    setEquationFlow(DNpr2);
    rVls,vars,nonliFuns,nonliVars=solve(DNpr2)
    somethingUpdated,fullDetermined=update!(DNpr2,rVls,vars)
  end
  #ideal gas solution 4738.739992314279
  println("result is:",4738.739992314279-DNpr1.h_Dep/1000.0+DNpr2.h_Dep/1000.0,"j/mol REF[2] p285 :3522")
  # REF[2] dh=3.522 (kj/mol)
end
function testPR()
  h_Dep::Array{Float64,1}=Array(Float64,0)
  v_Calc::Array{Float64,1}=Array(Float64,0)
  ###### verification: check thermodynamic prop of acetone Ref[1]:Table(2-185) for p=0.1Mpa #######
  ref_h=[47.730,49.643,54.255,59.166,64.436,70.066]*1000
  ref_s=[0.16988,0.17552,0.18783,0.19939,0.21049,0.22122]*1000
  ref_valuse=[25.930,28.002,32.563,36.923,41.200,45.437]
  ii=1
  for T in [328.84,350,400,450,500,550]
    DNpr=DANAPengRobinson()
    set(DNpr.pi)
    set(DNpr.R)
    set(DNpr.T,T)
    set(DNpr.P,0.1e6)
    # acetone
    tc,pc,af=getValueForCasNo("Criticals","67-64-1");
    set(DNpr.Tc,tc)
    set(DNpr.Pc,pc)
    set(DNpr.af,af)
    somethingUpdated=true
    fullDetermined=false
    nonliFuns::Array{Function,1}=Array(Function,0)
    nonliVars::Array{Set{String},1}=Array(Set{String},0)
    while (somethingUpdated && !fullDetermined)
      while (somethingUpdated && !fullDetermined)
        setEquationFlow(DNpr);
        rVls,vars,nonliFuns,nonliVars=solve(DNpr)
        println(DNpr.equationsFlow)
        somethingUpdated,fullDetermined=update!(DNpr,rVls,vars)
        return DNpr
      end
      if !fullDetermined
        i=1
        fullDetermined=true
        while (i<=length(nonliFuns))
          if length(nonliVars[i])==1
            result=Roots.fzero(nonliFuns[i],[0,typemax(Int64)])
            HelperEquation.setfield!(DNpr,nonliVars[i][1],result)
            somethingUpdated=true
          else
            fullDetermined=false
          end
          i=i+1
        end
      end
    end
    #println("T=",DNpr.T," v=",DNpr.v," ref_val=",ref_valuse[ii]) #Table[2.185]
    #println(" Dh=",DNpr.h_Dep," Ds=",DNpr.s_Dep) #Table[2.185]
    push!(h_Dep,get(DNpr.h_Dep))
    push!(v_Calc,get(DNpr.v))
    ii=ii+1
  end
  # ideal gas
  #println(h_Dep)
  #ideal_h=test_IdealGasEos.forAcetone()
  #res1=(ideal_h/1000.0+h_Dep/1000.0)
  #println("ideal gas dh (calculated)-> ",round((ideal_h-ideal_h[1])/1000.0)[2:end]," (j/mol)")
  #println("real gas dh (claculated)->  ",round(res1-res1[1])[2:end]," (j/mol)")
  #println("reference values->          ",round(ref_h-ref_h[1])[2:end]," (j/mol)")
end
=#
