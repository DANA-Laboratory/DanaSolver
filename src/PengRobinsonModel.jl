# REF[1] Engineering and Chemical Thermodynamics, 2nd Edition, Milo D. Koretsky
# REF[2] http://en.wikipedia.org/wiki/Departure_function
# REF[3] http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
# REF[4] Thermo-Hydro-Mechanical-Chemical Processes-in Fractured Porous Media
module PengRobinsonModel
  # Units J,Kmol,Kelvin,pascal
  using DanaTypes
  using EMLtypes
  import DanaTypes.setEquationFlow
  export DANAPengRobinson,setEquationFlow
  const R=8314.4621 #general gas constatnt "J/Kmol/Kelvin"
  const AVR_Tc=548.33512173913 #REF: Pkg.test("ThermodynamicsTable")
  const AVR_Pc=8.752551304347826e6
  const MAX_Tc=1113.0
  const MIN_Tc=5.2
  const MAX_Pc=6.162e8
  const MIN_Pc=227500.0
  type  DANAPengRobinson <: DanaModel
      DANAPengRobinson()=begin
        new(
          constant(Dict{Symbol,Any}(:Default=>pi)), 
          constant(Dict{Symbol,Any}(:Brief=>"general gas constatnt",:Default=>R,:Unit=>"J/Kmol/Kelvin")),
          "", 
          volume_mol(Dict{Symbol,Any}(:Upper=>50.0)),
          temperature(Dict{Symbol,Any}(:Upper=>1000.0)), 
          temperature(Dict{Symbol,Any}(:Brief=>"critical temperature")), 
          pressure(),
          pressure(Dict{Symbol,Any}(:Brief=>"critical pressure")),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(Dict{Symbol,Any}(:Lower=>eps(Float64),:Upper=>1/0.414,:Default=>0.1)), #B log(x>0)
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
          enth_mol(),
          entr_mol(),
          coefficient(Dict{Symbol,Any}(:Default=>0.457235*R^2*AVR_Tc^2/AVR_Pc,:Lower=>0.457235*R^2*MIN_Tc^2/MAX_Pc,:Upper=>0.457235*R^2*MAX_Tc^2/MIN_Pc)),
          [
            :(teta=acos(r/q^1.5)),
            :(Z1=-2*sqrt(q)*cos(teta/3)-beta/3),
            :(Z2=-2*sqrt(q)*cos((teta+2*pi)/3)-beta/3),
            :(Z3=-2*sqrt(q)*cos((teta+4*pi)/3)-beta/3),
            :(P=R*T/(v-b)-(a*(1+k*(1-sqrt(T/Tc)))^2)/(v*v+2*b*v-b*b)),
            :(Z=P*v/R/T),
            :(A=(a*(1+k*(1-sqrt(T/Tc)))^2)*P/R^2/T^2),#alpha=(1+k*(1-sqrt(T/Tc)))^2
            :(b=0.077796*R*Tc/Pc), #REF[3]
            :(B=b*P/R/T),
            :(beta=B-1),
            :(gama=A-3*B^2-2*B),
            :(delta=B^3+B^2-A*B),
            :(q=(beta*beta-3*gama)/9),
            :(r=(2*beta^3-9*beta*gama+27*delta)/54),
            :(s_Dep=R*(log(Z-B)-2.078*k*((1+k)/sqrt(T/Tc)-k)*log((Z+2.414*B)/(Z-0.414*B)))), #REF[2]
            :(h_Dep=R*T*(1-(v*((R*T)/(-b+v)-(a*(1+k*(1-sqrt(T/Tc)))^2)/(-b^2+2*b*v+v^2)))/(R*T)+(a*(1+k)*sqrt(T/Tc)*(-(k*T)+sqrt(T/Tc)*Tc+k*sqrt(T/Tc)*Tc)*(-log(-1+(b+v)/(sqrt(2)*b))+log(1+(b+v)/(sqrt(2)*b))))/(2*sqrt(2)*b*R*T^2))),
            :(a=0.457235*R^2*Tc^2/Pc), #REF[3]
            :(o=cbrt((r^2-q^3)^0.5+abs(r))),
            :(Z=-sign(r)*(o+q/o)-beta/3),
            :(k=0.37464+1.54226*af-0.26992*af^2), #REF[4] 2.87
            :(k=0.379642+1.48503*af-0.164423*af^2+0.016666*af^3a) #REF[4] 2.88
          ],Array(Expr,0)
        )
      end
      #paremeters
      pi::constant
      R::constant
      CASNO::String
      #variables
      v::volume_mol
      T::temperature
      Tc::temperature
      P::pressure
      Pc::pressure
      k::coefficient
      A::coefficient
      b::coefficient
      B::coefficient
      beta::coefficient
      gama::coefficient
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
      h_Dep::enth_mol
      s_Dep::entr_mol
      a::coefficient
      #equations
      equations::Array{Expr,1}
      equationsFlow::Array{Expr,1}
  end
  function setEquationFlow(this::DANAPengRobinson)
    if !(this.q.unset) && !(this.r.unset) && (get(this.q)^3-get(this.r)^2)<0
      this.equationsFlow=this.equations[5:19];
    else
      this.equationsFlow=this.equations[1:17];
    end
    if get(this.af)>0.49
      this.equationsFlow=[this.equationsFlow,this.equations[21]]
    else
      this.equationsFlow=[this.equationsFlow,this.equations[20]]    
    end

    if this.Z.unset && !(this.Z1.unset) && !(this.Z2.unset) && !(this.Z3.unset)
      set(this.Z,max(get(this.Z1),get(this.Z2),get(this.Z3)))
    end
  end
end
