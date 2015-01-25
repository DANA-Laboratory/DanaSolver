# REF[1] Engineering and Chemical Thermodynamics, 2nd Edition, Milo D. Koretsky
# REF[2] http://en.wikipedia.org/wiki/Departure_function
# REF[3] http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
module PengRobinsonModel
  # Units J,Kmol,Kelvin,pascal
  using DanaTypes
  using EMLtypes
  export DANAPengRobinson,setEquationFlow
  type  DANAPengRobinson <: DanaModel
      DANAPengRobinson()=begin
        new(
          constant(Dict(Symbol,Any)(:Default=>pi)),
          constant(Dict(Symbol,Any)(:Brief=>"general gas constatnt",:Default=>8314.4621,:Unit=>"J/Kmol/Kelvin")),
          "",
          volume_mol(),
          temperature(),
          temperature(Dict(Symbol,Any)(:Brief=>"critical temperature")),
          pressure(),
          pressure(Dict(Symbol,Any)(:Brief=>"critical pressure")),
          constant(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          coefficient(),
          constant(Dict(Symbol,Any)(:Brief=>"acentric factor")),
          coefficient(Dict(Symbol,Any)(:Brief=>"compressibility factor",:Lower=>eps(Float64),:Upper=>50.0)),
          coefficient(Dict(Symbol,Any)(:Brief=>"compressibility factor",:Lower=>eps(Float64),:Upper=>50.0)),
          coefficient(Dict(Symbol,Any)(:Brief=>"compressibility factor",:Lower=>eps(Float64),:Upper=>50.0)),
          coefficient(Dict(Symbol,Any)(:Brief=>"compressibility factor",:Lower=>eps(Float64),:Upper=>50.0)),
          coefficient(),
          enth_mol(),
          entr_mol(),
          enth_mol(),
          coefficient(Dict(Symbol,Any)(:Default=>0.457235*R^2*AVR_Tc^2/AVR_Pc,:Lower=>0.457235*R^2*MIN_Tc^2/MAX_Pc,:Upper=>0.457235*R^2*MAX_Tc^2/MIN_Pc)),
          [
            :(teta=acos(r/q^1.5)),
            :(Z1=-2*sqrt(q)*cos(teta/3)-beta/3),
            :(Z2=-2*sqrt(q)*cos((teta+2*pi)/3)-beta/3),
            :(Z3=-2*sqrt(q)*cos((teta+4*pi)/3)-beta/3),
            :(P=R*T/(v-b)-(a*(1+k*(1-sqrt(T/Tc)))^2)/(v*v+2*b*v-b*b)),
            :(Z=P*v/R/T),
            :(A=(a*(1+k*(1-sqrt(T/Tc)))^2)*P/R^2/T^2),
            :(b=0.077796*R*Tc/Pc), #REF[3]
            :(B=b*P/R/T),
            :(beta=B-1),
            :(gama=A-3*B^2-2*B),
            :(delta=B^3+B^2-A*B),
            :(q=(beta*beta-3*gama)/9),
            :(r=(2*beta^3-9*beta*gama+27*delta)/54),
            :(h_Dep=R*Tc*((T/Tc)*(Z-1)-2.078*(1+k)*sqrt((1+k*(1-sqrt(T/Tc)))^2)*log((Z+2.414*B)/(Z-0.414*B)))), #REF[2]
            :(s_Dep=R*(log(Z-B)-2.078*k*((1+k)/sqrt(T/Tc)-k)*log((Z+2.414*B)/(Z-0.414*B)))), #REF[2]
            :(h_Dep2=((-4*(b^3*R*T*Tc-2*b^2*R*T*Tc*v+a*(Tc-2*k*(-1+sqrt(T/Tc))*Tc+k^2*(T+Tc-2*sqrt(T/Tc)*Tc))*v^2-b*v*(a*(Tc-2*k*(-1+sqrt(T/Tc))*Tc+k^2*(T+Tc-2*sqrt(T/Tc)*Tc))+R*T*Tc*v)))/(Tc*(b-v)*(b^2-2*b*v-v^2))-(sqrt(2)*a*(1+k)*(-1+k*(-1+sqrt(T/Tc)))*log(-1+(b+v)/(sqrt(2)*b)))/b+(sqrt(2)*a*(1+k)*(-1+k*(-1+sqrt(T/Tc)))*log(1+(b+v)/(sqrt(2)*b)))/b)/4),
            :(a=0.457235*R^2*Tc^2/Pc), #REF[3]
            :(o=cbrt((r^2-q^3)^0.5+abs(r))),
            :(Z=-sign(r)*(o+q/o)-beta/3),
            :(k=0.37464+1.54226*af-0.26992*af^2),
            :(k=0.379642+1.48503*af-0.164423*af^2+0.016666*af^3)
            #:(PTv=-(((11431*k^2*Tc*v-11431*b*k^2*Tc)*R^2+(-25000*Pc*v^2-50000*b*Pc*v-25000*b^2*Pc)*R)*sqrt(T)+sqrt(Tc)*((-11431*k^2-11431*k)*Tc*v+(11431*b*k^2+11431*b*k)*Tc)*R^2)/((25000*Pc*v^3+25000*b*Pc*v^2-75000*b^2*Pc*v+25000*b^3*Pc)*sqrt(T))),
            #s calculation REF[1] EQ(5.31)
            #:(PTvIv=R*(sqrt((2)*(50000*b*Pc*sqrt(T)-11431*k*R*(-(k*sqrt(T))+sqrt(Tc)+k*sqrt(Tc))*Tc)*atanh((b+v)/(sqrt(2)*b))+25000*b*Pc*sqrt(T)*(4*log(-b+v)-log(-b^2+2*b*v+v^2))))/(50000*b*Pc*sqrt(T))),
            #u calculation REF[1] EQ(5.36)
            #:(TmuPTvmiPIV=-(P*v)+(R*sqrt(T)*(50000*b*Pc*Sqrt(T)+11431*k^2*R*Sqrt(T)*Tc-11431*k*R*Tc^(3/2)-11431*k^2*R*Tc^(3/2))*atanh((b+v)/(sqrt(2)*b)))/(25000*sqrt(2)*b*Pc)+2*R*T*log(-b+v)-(R*T*log(-b^2+2*b*v+v^2))/2),
            #:(PvT=-(((-11431*k^2*Tc*v^3+11431*b*k^2*Tc*v^2+11431*b^2*k^2*Tc*v-11431*b^3*k^2*Tc)*R^2+(12500*Pc*v^4+50000*b*Pc*v^3+100000*b^2*Pc*v^2-62500*b^4*Pc)*R)*T+sqrt(Tc)*((22862*k^2+22862*k)*Tc*v^3+(-22862*b*k^2-22862*b*k)*Tc*v^2+(-22862*b^2*k^2-22862*b^2*k)*Tc*v+(22862*b^3*k^2+22862*b^3*k)*Tc)*R^2*sqrt(T)+((-11431*k^2-22862*k-11431)*Tc^2*v^3+(11431*b*k^2+22862*b*k+11431*b)*Tc^2*v^2+(11431*b^2*k^2+22862*b^2*k+11431*b^2)*Tc^2*v+(-11431*b^3*k^2-22862*b^3*k-11431*b^3)*Tc^2)*R^2)/(12500*Pc*v^6+25000*b*Pc*v^5-62500*b^2*Pc*v^4-50000*b^3*Pc*v^3+137500*b^4*Pc*v^2-75000*b^5*Pc*v+12500*b^6*Pc)),
            #:(PvTIT=(R*(b+v)*(11431*b^2*(1+k)^2*R*Tc^2*T-22862*b*(1+k)^2*R*Tc^2*v*T+11431*(1+k)^2*R*Tc^2*v^2*T-(45724*b^2*k*(1+k)*R*Tc^(3/2)*T^(3/2))/3+(91448*b*k*(1+k)*R*Tc^(3/2)*v*T^(3/2))/3-(45724*k*(1+k)*R*Tc^(3/2)*v^2*T^(3/2))/3+31250*b^3*Pc*T^2+(b^2*(11431*k^2*R*Tc-62500*Pc*v)*T^2)/2+(v^2*(11431*k^2*R*Tc-12500*Pc*v)*T^2)/2-b*v*(11431*k^2*R*Tc+18750*Pc*v)*T^2))/(12500*Pc*(b^3-3*b^2*v+b*v^2+v^3)^2)),
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
      k::constant
      A::coefficient
      b::coefficient
      B::coefficient
      beta::coefficient
      gama::coefficient
      delta::coefficient
      q::coefficient
      r::coefficient
      teta::coefficient
      af::constant
      Z1::coefficient
      Z2::coefficient
      Z3::coefficient
      Z::coefficient
      o::coefficient
      h_Dep::enth_mol
      s_Dep::entr_mol
      h_Dep2::enth_mol
      a::coefficient
      #h_Dep_aletr::Float64
      #equations
      equations::Array{Expr,1}
      equationsFlow::Array{Expr,1}
  end
  function setEquationFlow(this::DANAPengRobinson)
    if !(this.q.unset) && !(this.r.unset) && (get(this.q)^3-get(this.r)^2)<0
      this.equationsFlow=this.equations[5:20];
    else
      this.equationsFlow=this.equations[1:18];
    end
    if get(this.af)>0.49
      this.equationsFlow=[this.equationsFlow,this.equations[22]]
    else
      this.equationsFlow=[this.equationsFlow,this.equations[21]]    
    end

    if this.Z.unset && !(this.Z1.unset) && !(this.Z2.unset) && !(this.Z3.unset)
      set(this.Z,max(get(this.Z1),get(this.Z2),get(this.Z3)))
    end
  end
end
