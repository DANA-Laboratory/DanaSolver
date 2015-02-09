# REF[1] Engineering and Chemical Thermodynamics, 2nd Edition, Milo D. Koretsky
# REF[2] http://en.wikipedia.org/wiki/Departure_function
# REF[3] http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
# REF[4] Thermo-Hydro-Mechanical-Chemical Processes-in Fractured Porous Media
module PengRobinsonDLModel
  # Units J,Kmol,Kelvin,pascal
  using DanaTypes
  using EMLtypes
  import DanaTypes.setEquationFlow
  export DANAPengRobinsonDL,setEquationFlow
  const R=8314.4621 #general gas constatnt "J/Kmol/Kelvin"
  const AVR_Tc=548.33512173913 #REF: Pkg.test("ThermodynamicsTable")
  const AVR_Pc=8.752551304347826e6
  const MAX_Zc=0.428
  const MIN_Zc=0.117
  const MAX_Tc=1113.0
  const MIN_Tc=5.2
  const MAX_Pc=6.162e8
  const MIN_Pc=227500.0
  type  DANAPengRobinsonDL <: DanaModel
      DANAPengRobinsonDL()=begin
        new(
          constant(Dict{Symbol,Any}(:Default=>pi)), 
          constant(Dict{Symbol,Any}(:Brief=>"general gas constatnt",:Default=>R,:Unit=>"J/Kmol/Kelvin")),
          constant(Dict{Symbol,Any}(:Brief=>"critical compressiblity for pr",:Default=>0.31115542517055555)),
          "", 
          #parameters
          coefficient(Dict{Symbol,Any}(:Brief=>"critical compressibility",:Lower=>eps(Float64),:Upper=>1.0,:Default=>0.8)),
          temperature(Dict{Symbol,Any}(:Brief=>"critical temperature")), 
          pressure(Dict{Symbol,Any}(:Brief=>"critical pressure")),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(Dict{Symbol,Any}(:Brief=>"acentric factor")),
          coefficient(Dict{Symbol,Any}(:Brief=>"compressibility factor",:Lower=>eps(Float64),:Upper=>1.0,:Default=>0.8)),
          coefficient(Dict{Symbol,Any}(:Brief=>"compressibility factor",:Lower=>eps(Float64),:Upper=>1.0,:Default=>0.8)),
          coefficient(Dict{Symbol,Any}(:Brief=>"compressibility factor",:Lower=>eps(Float64),:Upper=>1.0,:Default=>0.8)),
          coefficient(Dict{Symbol,Any}(:Brief=>"compressibility factor",:Lower=>eps(Float64),:Upper=>1.0,:Default=>0.8)),
          coefficient(),
          coefficient(),
          #variables
          coefficient(Dict{Symbol,Any}(:Brief=>"volume dimensionless")),
          coefficient(Dict{Symbol,Any}(:Brief=>"temprature dimensionless")),
          coefficient(Dict{Symbol,Any}(:Brief=>"pressure dimensionless")),
          coefficient(Dict{Symbol,Any}(:Brief=>"enthalpy dimensionless")),
          coefficient(Dict{Symbol,Any}(:Brief=>"entropy dimensionless")),
          coefficient(Dict{Symbol,Any}(:Brief=>"gibbs dimensionless")),
          [
            :(teta=acos(r/q^1.5)),
            :(Z1=-2*sqrt(q)*cos(teta/3)-beta/3),
            :(Z2=-2*sqrt(q)*cos((teta+2*pi)/3)-beta/3),
            :(Z3=-2*sqrt(q)*cos((teta+4*pi)/3)-beta/3),
            :(Pr=Tr/(vr-br)/Zc-(ar*alpha)/(vr^2+2*br*vr-br^2)),
            :(Z=Pr*vr*Zc/Tr),
            :(A=(ar*alpha)*Pr*Zc^2/Tr^2),
            :(alpha=(1+k*(1-sqrt(Tr)))^2),
            :(br=0.077796/Zc),#br=b/vc #REF[3]
            :(B=br*Zc*Pr/Tr),
            :(beta=B-1),
            :(gama=A-3*B^2-2*B),
            :(delta=B^3+B^2-A*B),
            :(q=(beta*beta-3*gama)/9),
            :(r=(2*beta^3-9*beta*gama+27*delta)/54),
            #0.457235/0.077796/2/sqrt(2)=2.0779601078193677~2.07796
            :(s_Dep=2.0779601078193677*k*((1+k)/sqrt(Tr)-k)*log(((1+sqrt(2))*br+vr)/((1-sqrt(2))*br+vr))-log(Z-B)),
            :(h_Dep=1-Z+(2.0779601078193677*sqrt(alpha)*(1+k)*log(((1+sqrt(2))*br+vr)/((1-sqrt(2))*br+vr)))),
            :(g_Dep=1-Z+(alpha*ar*Zc*log(((1+sqrt(2))*br+vr)/((1-sqrt(2))*br+vr)))/(sqrt(2)*br*Tr)+log(((-br+vr)^Zc*Z)/vr)),
            :(ar=0.457235/Zc^2), #REF[3]
            :(o=cbrt((r^2-q^3)^0.5+abs(r))),
            :(Z=-sign(r)*(o+q/o)-beta/3),
            :(k=0.37464+1.54226*af-0.26992*af^2), #REF[4] 2.87
            :(k=0.379642+1.48503*af-0.164423*af^2+0.016666*af^3) #REF[4] 2.88
          ],
          Array(Expr,0),
          [:pi,:R,:name],
          [:Zc,:Tc,:Pc,:k,:beta,:gama,:alpha,:A,:B,:br,:delta,:q,:r,:teta,:af,:Z1,:Z2,:Z3,:Z,:o,:ar],
          [:vr,:Tr,:Pr,:h_Dep,:s_Dep,:g_Dep]
        )
      end
      #consts
      pi::constant
      R::constant
      Zc::constant
      name::String
      #parameters
      Tc::temperature
      Pc::pressure
      k::coefficient
      beta::coefficient
      gama::coefficient
      alpha::coefficient
      A::coefficient
      B::coefficient
      br::coefficient
      delta::coefficient
      q::coefficient
      r::coefficient
      teta::coefficient
      af::coefficient
      Z1::coefficient
      Z2::coefficient
      Z3::coefficient
      Z::coefficient
      o::coefficient
      ar::coefficient
      #variables
      vr::coefficient
      Tr::coefficient
      Pr::coefficient
      h_Dep::coefficient
      s_Dep::coefficient
      g_Dep::coefficient
      #equations
      equations::Array{Expr,1}
      equationsFlow::Array{Expr,1}
      constants::Array{Symbol,1}
      parameters::Array{Symbol,1}
      variables::Array{Symbol,1}
  end
  function setEquationFlow(this::DANAPengRobinsonDL)
    if !(this.q.unset) && !(this.r.unset) && (get(this.q)^3-get(this.r)^2)<0
      this.equationsFlow=this.equations[5:21];
    else
      this.equationsFlow=this.equations[1:19];
    end
    if get(this.af)>0.49
      this.equationsFlow=[this.equationsFlow,this.equations[22]]
    else
      this.equationsFlow=[this.equationsFlow,this.equations[23]]    
    end

    if this.Z.unset && !(this.Z1.unset) && !(this.Z2.unset) && !(this.Z3.unset)
      set(this.Z,max(get(this.Z1),get(this.Z2),get(this.Z3)))
    end
  end
end
