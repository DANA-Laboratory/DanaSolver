
@test_approx_eq Solver.callfzero(x->x^3-1,0.0,[0.0,100.0])[1] 1.0

@time r=Solver.callfzero(x->x^2-1,200.0,[0.0,500.0])
println("attempts=$r")
@test_approx_eq Solver.callfzero(x->x^2-1,0.0,[0.0,100.0])[1] 1.0

@time r=Solver.callfzero(x->x^3+log(2.0-x)-1,50.0,[0.0,100.0],1000)
@test_approx_eq r[1] 1.0
println("attempts=$r")

@test_throws ErrorException Solver.callfzero(x->x^2+log(2.0-x)-1,-50.0,[-100.0,-10.0],100)

@test_throws DomainError Solver.callfzero(x->x^2+log(2.0-x)-1,500.0,[100.0,10000.0],100)

@time r=Solver.callfzero(y->(y)^3+log(y-1),50.0,[0.0,100.0],1000)
println("(y)^3+log(y-1) ->y=",r)

println("********test opt********")
f(x,y)=begin
        println("in f $x , $y");
        return x^3+log(y)
       end
g(x,y)=begin
        println("in g $x , $y");
        return log(x-y)
       end
indG=Array((Vector{Int}),0)
push!(indG,[1,2],[1,2])
# cant solve because it sends y>x and log(x-y) fails 
@test_throws DomainError Solver.callmin([f,g],indG,[0.0,0.0],[100.0,100.0],[10.0,1.0])

println("********test fzero(f(x,fzero(g(x,y)))) instead of opt ********")
@time r=Solver.callfzero(x->f(x,Solver.callfzero(y->g(x,y),50.0,[0.0,100.0],1000)[1]),50.0,[0.0,100.0],1000)
println("use 2 fzeroz ->",r)

function callffzero(x,g::Function,gde,gbracket,f::Function,fde,fbracket)
  gde=Solver.callfzero(g,gde,gbracket,1000)[1]
  return f(x,fde,fbracket,1000)
end
#=
println("********test fzero(g(x,fzero(f(x,y)))) instead of opt ********")
@time r=Solver.callfzero(x->g(x,Solver.callfzero(y->f(x,y),50.0,[0.0,100.0],1000)[1]),50.0,[0.0,100.0],1000)
println("use 2 fzeroz ->",r)
=#