type Equation
  Equation(ex::Expr)=new(ex,Set{String}(),Set{String}(),Array(Float64,0))
  ex::Expr
  termall::Set{String} #terms of a equation
  termnonli::Set{String} #non-linear terms of a equation
  factliinnonli::Vector{Float64} #linear terms in nolinear equations
end
type Equations
  Equations(exs::Array{Equation,1})=new(exs,Array(Float64,0,0),Array(String,0),Array(Int,0))
  exs::Vector{Equation}
  facts::Matrix{Float64} #matrix of factors of linear equations
  terms::Vector{String} #array of symbols (all term of eqautions)
  indexnonliexs::Vector{Int} #index of nonlinear equations
end
