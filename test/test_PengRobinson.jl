# REF[1] Chemical Process Design and Integration By Robin Smith
# REF[2] Engineering and Chemical Thermodynamics By Milo D. Koretsky
# REF[3] Perry HandBook
reload ("HelperEquation.jl")
reload ("Tables.jl")
reload ("PengRobinson.jl")
reload ("Calculus.jl")
#reload("PengRobinson.jl/test/test_PengRobinson.jl");test_PengRobinson.testPR()
#reload("PengRobinson.jl/test/test_PengRobinson.jl");test_PengRobinson.testDeparture()
#reload("PengRobinson.jl/test/test_PengRobinson.jl");test_PengRobinson.testVariousKnowns()
#reload("PengRobinson.jl/test/test_PengRobinson.jl");test_PengRobinson.testMoreThanOneNonLinear()
module test_PengRobinson
  using HelperEquation
  using PengRobinson
  using Tables
  using Roots
	using NLopt
	using DanaTypes
  reload ("IdealGasEos.jl/test/test_IdealGasEos.jl")
  using test_IdealGasEos
	#generate indexes for all possible system of _noe equations[APSOE]. 
  #where equations is a selection from _minIndex to _maxIndex of a list of equations
	function getAPSOE(_minIndex::Int,_maxIndex::Int,_noe::Int)
		if 1<_noe
			jj::Vector=Vector[]
			for k in [_minIndex+1:_maxIndex] 
				j=getAPSOE(k,_maxIndex,_noe-1)
				jj=append!(jj,[push!(e,k-1) for e in j])
			end
			return jj
		else
			return [[i] for i in _minIndex:_maxIndex]
		end
	end
	#*********************
	function testMoreThanOneNonLinear()
		nonliFuns::Array{Function,1}=Array(Function,0)
		nonliVars::Array{Set{String},1}=Array(Set{String},0)
		for k in [1:5]
			cNo="75-07-0"
			PR=DANAPengRobinson()
			PR.Tc,PR.Pc,PR.af=getValueForCasNo("Criticals",cNo) 
			PR.h_Dep,PR.h_Dep2,PR.P,PR.T,PR.v=[(NaN,NaN,PR.Pc,PR.Tc,NaN),(NaN,NaN,NaN,PR.Tc,0.22435958119569624),(NaN,NaN,PR.Pc,NaN,0.22435958119569624),(NaN,-1.08307222504879e7,PR.Pc,NaN,NaN),(-1.08307222504879e7,NaN,PR.Pc,NaN,NaN)][k]
			somthingUpdated=true
			fullDetermined=false
			while (somthingUpdated && !fullDetermined)
				while (somthingUpdated && !fullDetermined)
					while (somthingUpdated && !fullDetermined)
						setEquationFlow(PR);
						#linear solver
						rVls,vars,nonliFuns,nonliVars=solve(PR)
						somthingUpdated,fullDetermined=update!(PR,rVls,vars)
						println("Linear Solver.")
					end
					if !fullDetermined
						somthingUpdated=false
						i=1
						fullDetermined=true
						#search for non-linear equations with only one unknown
						numberOfEquations=length(nonliFuns)
						while (i<=numberOfEquations && !somthingUpdated)
							if length(nonliVars[i])==1
								result=Roots.fzero(nonliFuns[i],[0,typemax(Int64)])
								HelperEquation.setfield(PR,[nonliVars[i]...][1],result)
								somthingUpdated=true
								println("non-Linear Solver, solves one equation.")
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
					if (!somthingUpdated)
						i=2
						while (i<=numberOfEquations && !somthingUpdated)
							eqIndexes=getAPSOE(1,numberOfEquations,i)
							for eqIndex in eqIndexes
								varGroup=getindex(nonliVars,eqIndex)
								allVars=union(varGroup...)
								if length(allVars) == i
									eqGroup=getindex(nonliFuns,eqIndex)
									println("eqGroup=",eqIndex," for vars:",allVars)
									indxGroup=map(x->indexin([x...],[allVars...]),varGroup)
									opt = Opt(:GN_DIRECT_L, i)
									lower_bounds!(opt, [1.0e-3, 1.0])
									upper_bounds!(opt, [10.0,2500])
									stopval!(opt, 1.0e-12)
									maxtime!(opt, 1.0*i)
									#ftol_abs!(opt, 1.0e-19)
									#ftol_rel!(opt, 1.0e-18)
									min_objective!(opt, (y,gradient)->mapreduce(x->(apply(eqGroup[x],getindex(y,indxGroup[x])))^2,+,[1:i]))
									(minf,minx,ret)=optimize(opt,ones(Float64,i))
									println("got $minf at $minx (returned $ret)")
									if "$ret"=="STOPVAL_REACHED"
										for j in [1:i]
											HelperEquation.setfield(PR,[nonliVars[eqIndex[1]]...][j],minx[j])
										end
										somthingUpdated=true
										if i==numberOfEquations
											fullDetermined=true
										end
										break
									end
								end
							end
							i+=1
						end
					end
				end
			end
			if fullDetermined
				println("Solution Done! v=",round(PR.v,7)," T=",round(PR.T,7)," P=",round(PR.P,7)," h=",round(PR.h_Dep2,7))
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
		somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(PR);
      rVls,vars,nonliFuns,nonliVars=solve(PR)
      somthingUpdated,fullDetermined=update!(PR,rVls,vars)
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
		somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(PR);
      rVls,vars,nonliFuns,nonliVars=solve(PR)
      somthingUpdated,fullDetermined=update!(PR,rVls,vars)
			if fullDetermined
				println("solved for P! PR.P=",PR.P)
			end
		end
		#v & P
		PR=DANAPengRobinson()
		PR.Tc,PR.Pc,PR.af=getValueForCasNo("Criticals",cNo)
		PR.v=v
		PR.P=PR.Pc
		somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(PR);
      rVls,vars,nonliFuns,nonliVars=solve(PR)
      somthingUpdated,fullDetermined=update!(PR,rVls,vars)
			if !fullDetermined
				i=1
				fullDetermined=true
				while (i<=length(nonliFuns))
					if length(nonliVars[i])==1
						result=Roots.fzero(nonliFuns[i],[0,typemax(Int64)])
						HelperEquation.setfield(PR,nonliVars[i][1],result)
						somthingUpdated=true
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
    res=optimize(optFunctionHP, [PR.T/20,PR.h_Dep2/20]);
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
    h_Dep2=x[2];
    ret1=h_Dep2-((-4*(b^3*R*T*Tc-2*b^2*R*T*Tc*v+d*(Tc-2*k*(-1+sqrt(T/Tc))*Tc+k^2*(T+Tc-2*sqrt(T/Tc)*Tc))*v^2-b*v*(d*(Tc-2*k*(-1+sqrt(T/Tc))*Tc+k^2*(T+Tc-2*sqrt(T/Tc)*Tc))+R*T*Tc*v)))/(Tc*(b-v)*(b^2-2*b*v-v^2))-(sqrt(2)*d*(1+k)*(-1+k*(-1+sqrt(T/Tc)))*log(-1+(b+v)/(sqrt(2)*b)))/b+(sqrt(2)*d*(1+k)*(-1+k*(-1+sqrt(T/Tc)))*log(1+(b+v)/(sqrt(2)*b)))/b)/4;
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
		somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(DNpr1);
      rVls,vars,nonliFuns,nonliVars=solve(DNpr1)
      somthingUpdated,fullDetermined=update!(DNpr1,rVls,vars)
    end
    somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(DNpr2);
      rVls,vars,nonliFuns,nonliVars=solve(DNpr2)
      somthingUpdated,fullDetermined=update!(DNpr2,rVls,vars)
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
			somthingUpdated=true
			fullDetermined=false
			nonliFuns::Array{Function,1}=Array(Function,0)
			nonliVars::Array{Set{String},1}=Array(Set{String},0)
			while (somthingUpdated && !fullDetermined)
				while (somthingUpdated && !fullDetermined)
          setEquationFlow(DNpr);
					rVls,vars,nonliFuns,nonliVars=solve(DNpr)
					println(DNpr.equationsFlow)
					somthingUpdated,fullDetermined=update!(DNpr,rVls,vars)					
					return DNpr
				end
				if !fullDetermined
					i=1
					fullDetermined=true
					while (i<=length(nonliFuns))
						if length(nonliVars[i])==1
							result=Roots.fzero(nonliFuns[i],[0,typemax(Int64)])
							HelperEquation.setfield(DNpr,nonliVars[i][1],result)
							somthingUpdated=true
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
end