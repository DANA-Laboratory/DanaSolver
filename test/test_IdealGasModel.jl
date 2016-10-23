# REF[2] Engineering and Chemical Thermodynamics By Milo D. Koretsky
# REF[2] Example 5.4
function forButane()
  println("**** test for butane ****")
  dh=0.0;h1=0.0;
  for T in [80+273.15,120+273.15]
    DNIdel=DANAIdealGasEos()
    DNIdel.name="Butane" #CASNO = getvalueforname("Profile","Butane")[2]
    DNIdel.usePolynomialEstimationOfCp=false
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpHyper",DNIdel.name)
    DNIdel.T=T
    setEquationFlow(DNIdel)
    somethingUpdated,fullDetermined=Solver.slstsubnfd!(DNIdel)
    dh=DNIdel.h-h1
    h1=DNIdel.h
  end
  println("result is=",dh/1000,"j/mol ref[2] page285:",4696)
end
function CompPolyVS_Hyper()
  println("**** butane CompPolyVS_Hyper ****")
  h_poly::Array{Float64,1}=Array(Float64,0)
  h_hyper::Array{Float64,1}=Array(Float64,0)
  for T in [298.15,328.84,350,400,450,500,550]
    DNIdelPoly=DANAIdealGasEos()
    DNIdelHyper=DANAIdealGasEos()
	  DNIdelPoly.T=T
		DNIdelHyper.T=T
    DNIdelPoly.name="Butane"
    DNIdelHyper.name="Butane"
    DNIdelPoly.usePolynomialEstimationOfCp=true
    DNIdelHyper.usePolynomialEstimationOfCp=false
    DNIdelPoly.C1,DNIdelPoly.C2,DNIdelPoly.C3,DNIdelPoly.C4,DNIdelPoly.C5 = getvalueforname("CpPoly","Butane")
    DNIdelHyper.C1,DNIdelHyper.C2,DNIdelHyper.C3,DNIdelHyper.C4,DNIdelHyper.C5 = getvalueforname("CpHyper","Butane")
    setEquationFlow(DNIdelPoly)
    setEquationFlow(DNIdelHyper)
    somethingUpdated,fullDetermined=Solver.slstsubnfd!(DNIdelPoly)
    somethingUpdated,fullDetermined=Solver.slstsubnfd!(DNIdelHyper)
    push!(h_poly,DNIdelPoly.h)
    push!(h_hyper,DNIdelHyper.h)
  end
  println("poly ->",h_poly-h_poly[1])
  println("hyper->",h_hyper-h_hyper[1])
end
function forwater()
  println("**** for water ****")
  i=1;
  h_Dep::Array{Float64,1}=Array(Float64,0)
  cp_Dep::Array{Float64,1}=Array(Float64,0)
	tem = [298.15,200,240,260,280,300,320,340,360,380,400,420,440,460,480,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3500,4000,4500,5000]
  for T in tem
	  DNIdel=DANAIdealGasEos()
		DNIdel.T=T
	  DNIdel.name="Water" #7732-18-5
		DNIdel.usePolynomialEstimationOfCp=false
		DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpHyper",DNIdel.name)
		setEquationFlow(DNIdel)
    somethingUpdated,fullDetermined=Solver.slstsubnfd!(DNIdel)
    push!(h_Dep,DNIdel.h)
    push!(cp_Dep,DNIdel.Cp)
		i+=1;
	end
  # println(" Cp=",cp_Dep);
  ref_h=[0,-3282,-1948,-1279,-609,62,735,1410,2088,2769,3452,4139,4829,5523,6222,6925,8699,10501,12321,14192,16082,18002,19954,21938,23954,26000,30191,34506,38942,43493,48151,52908,57758,62693,67706,72790,77941,83153,88421,93741,99108,104520,109973,115464,120990,126549,154768,183552,212764,242313]
  DhDep=(h_Dep-h_Dep[1])/1000
  i=1
  for t in tem
    println("at $t -> calculated enthalpy=",DhDep[i]," Ref[1] Table2-180=",ref_h[i]," %diff=",(DhDep[i]-ref_h[i])/ref_h[i]*100)
    i+=1
  end
end
function forCarbonmonoxide()
  println("****  for Carbonmonoxide ****")
  ###### verification: check monoxide Enthalpies with Ref[1]:Table(2-180) #######
  DNIdel=DANAIdealGasEos()
  DNIdel.T=298.15
  # monoxide
  DNIdel.name="Carbon monoxide" #"630-08-0"
  DNIdel.usePolynomialEstimationOfCp=false
  DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpHyper","Carbon monoxide")
  setEquationFlow(DNIdel)
  somethingUpdated,fullDetermined=Solver.slstsubnfd!(DNIdel)
  hst=DNIdel.h;
	ust=DNIdel.u;
	refDelta_h=[-2858,-1692,-1110,-529,54,638,1221,1805,2389,2975,3563,4153,4643,5335,5931,7428,8942,10477,12023,13592,15177,16781,18401,20031,21690,25035,28430,31868,35343,38850,42385,45945,49526,53126,56744,60376,64021,67683,71324,74985,78673,82369,86074,89786,93504,112185,130989,149895,168890];
	i=1;
	for T in [200,240,260,280,300,320,340,360,380,400,420,440,460,480,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3500,4000,4500,5000]	
		DNIdel=DANAIdealGasEos()
		DNIdel.T=T
		# monoxide
    DNIdel.name="Carbon monoxide"#"630-08-0"
		DNIdel.usePolynomialEstimationOfCp=false
		DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getvalueforname("CpHyper","Carbon monoxide")
		setEquationFlow(DNIdel)
    somethingUpdated,fullDetermined=Solver.slstsubnfd!(DNIdel)
		println("T=",T," Dh=",(DNIdel.h-hst)/1000," ref value=",refDelta_h[i]," %diff=",((DNIdel.h-hst)/1000-refDelta_h[i])/refDelta_h[i]*100);
		i+=1;
	end
end
