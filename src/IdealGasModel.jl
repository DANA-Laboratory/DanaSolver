# REF[1] Engineering and Chemical Thermodynamics, 2nd Edition, Milo D. Koretsky
# REF[2] Perry 8ed.
module IdealGasModel
  #Units J,Kmol,Kelvin,pascal
  using DanaTypes
  import DanaTypes.setEquationFlow
  export DANAIdealGasEos,setEquationFlow
  type  DANAIdealGasEos <: DanaModel
      DANAIdealGasEos()=begin
        new(8314.4621,"",true,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
          [
            :(P*v=R*T),#pascal
            :(Cp=C1+C2*T+C3*T^2+C4*T^3+C5*T^4),#Poly Cp in J/[Kmol*K] REF[2] TABLE 2-155
            :(Cp=C1+C2*((C3/T)/sinh(C3/T))^2+C4*((C5/T)/cosh(C5/T))^2),#Hyper Cp in J/[Kmol*K] REF[2] TABLE 2-156
            :(ICpOnTDT=C2*T+(C3*T^2)/2+(C4*T^3)/3+(C5*T^4)/4+C1*log(T)),#Integral of Cp/T Poly
            :(ICpOnTDT=(C2*C3*coth(C3/T)+C1*T*log(T)+C4*T*log(cosh(C5/T))-C2*T*log(sinh(C3/T))-C4*C5*tanh(C5/T))/T),#Integral of Cp/T Hyper
            :(Cv=Cp-R),#Cv Def
            :(ICpDT=C1*T+1/2*C2*T^2+1/3*C3*T^3+1/4*C4*T^4+1/5*C5*T^5),#Integ of Cp Poly
            :(ICpDT=C1*T+C2*C3*coth(C3/T)-C4*C5*tanh(C5/T)),#Integ of Cp Hyper
            :(u=ICpDT-R*T), #Internal energy in J/Kmol
            :(h=ICpDT), #Enthalpy in J/kmol
            :(s=ICpOnTDT-R*log(P)), #Entropy in J/[Kmol*K] ,REF[1] eq(3.22)
            :(f=u-T*s), #Helmholtz free energy
            :(g=h-T*s) #Gibbs free energy
          ],Array(Expr,0)
        )
      end
      #paremeters
      R::Float64
      name::String
      usePolynomialEstimationOfCp::Bool
      C1::Float64
      C2::Float64
      C3::Float64
      C4::Float64
      C5::Float64
      #variables
      v::Float64
      T::Float64
      P::Float64
      Cp::Float64
      Cv::Float64
      u::Float64
      h::Float64
      s::Float64
      g::Float64
      f::Float64
      ICpOnTDT::Float64
      ICpDT::Float64
      #equations
      equations::Array{Expr,1}
      equationsFlow::Array{Expr,1}
  end
  function setEquationFlow(this::DANAIdealGasEos)
    if this.usePolynomialEstimationOfCp
      this.equationsFlow=[this.equations[1],this.equations[2],this.equations[4],this.equations[6],this.equations[7],this.equations[9],this.equations[10],this.equations[11],this.equations[12],this.equations[13]]
    else
      this.equationsFlow=[this.equations[1],this.equations[3],this.equations[5],this.equations[6],this.equations[8],this.equations[9],this.equations[10],this.equations[11],this.equations[12],this.equations[13]]
    end
  end
end
