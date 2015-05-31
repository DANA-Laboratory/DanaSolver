# REF[2] Perry 8ed.
module LiquidsModel
using DanaTypes
import DanaTypes.setEquationFlow
export DANALiquids,setEquationFlow
type  DANALiquids <: DanaModel
  DANALiquids()=begin
    new(
    [
    :(density=Cd1/(Cd2^(1+(1-T/Cd3)^Cd4))),
    :(density=Cd1+Cd2*T+Cd3*T^2+Cd4*T^3),
    :(density=17.863+58.606t^0.35−95.396t^(2/3)+213.89t−141.26t^(4/3)),
    :(t=1−T/647.096),
    :(vp=exp(Cvp1+Cvp2/T+Cvp3*ln(T)+Cvp4*T^Cvp5)),
    :(cp=Ccp1+Ccp2*T+Ccp3*T^2+Ccp4*T^3+Ccp4*T^4),
    :(cp=(Ccp1^2)/tt+Ccp2−2*Ccp1*Ccp3*tt−Ccp1*Ccp4*tt^2−(Ccp3^2*tt^3)/3−(Ccp3*Ccp4*tt^4)/2−(Ccp4*tt^5)/5),
    :(tt=1-T/Tc)
    ],Array(Expr,0)
    )
  end
  #parameters
  name::String
  Cd1::Float64
  Cd2::Float64
  Cd3::Float64
  Cd4::Float64
  Tmind::Float64
  Tmaxd::Float64
  Cvp1::Float64
  Cvp2::Float64
  Cvp3::Float64
  Cvp4::Float64
  Cvp5::Float64
  Tminvp::Float64
  Tmaxvp::Float64
  Ccp1::Float64
  Ccp2::Float64
  Ccp3::Float64
  Ccp4::Float64
  Ccp5::Float64
  Tmincp::Float64
  Tmaxcp::Float64
  Tc::Float64
  #variables
  density::Float64
  t::Float64
  vp::Float64
  cp::Float64
  tt::Float64
  #equations
  equations::Array{Expr,1}
  equationsFlow::Array{Expr,1}
end
function setEquationFlow(this::DANALiquids)
  if !(name in ["o-Terphenyl [note: limited range]","Water [note: limited range]","Water"])
    this.equationsFlow=[this.equations[1]]
  else
    if(name=="Water")
      this.equationsFlow=[this.equations[3:4]]
    else
      this.equationsFlow=[this.equations[2]]
    end
  end
  if !(name in ["Ammonia [note: use Eq.(2)]","1,2-Butanediol [note: use Eq(2]","1,3-Butanediol [note: use Eq(2)]","Carbon monoxide [note: use Eq.(2)]","1,1-Difluoroethane [note: use Eq.(2)]","Ethane [note: use Eq.(2)]","Heptane [note: use Eq.(2)]","Hydrogen [note: use Eq.(2)]","Hydrogen sulfide [note: use Eq.(2)]","Methane [note: use Eq.(2)]","Propane [note: use Eq.(2)]"])
    this.equationsFlow=[this.equationsFlow,this.equations[6]]
  else
    this.equationsFlow=[this.equationsFlow,this.equations[7:8]]
  end
  this.equationsFlow=[this.equationsFlow,this.equations[5]]
end
end
