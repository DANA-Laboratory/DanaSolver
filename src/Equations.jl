type Equations
  eq::Array{Expr,1} #array of equations
  facts::Array{Float64,2} #matrix of factors of linear equations
  terms::Array{String,1} #array of symbols (all term of eqautions)
  indexnonlieq::Array{Int,1} #index of nonlinear equations
  termsli::Array{Set{String}i,1} #array of sets each expresses linear terms of a equation
  termsnonli::Array{Set{String}i,1} #array of sets each expresses non-linear terms of a equation
  factsliinnonli::Array{Array{Float64,1},1} #factors of linear terms in nolinear equations
end
