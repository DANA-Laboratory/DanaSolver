output=false

@test_approx_eq Solver.callfzero(x->x^3-1,0.0,[0.0,100.0])[1] 1.0

@time r=Solver.callfzero(x->x^2-1,200.0,[0.0,500.0])
@test r == (1.0,1)

@time r=Solver.callfzero(x->x^3+log(2.0-x)-1,50.0,[0.0,100.0],1000)
@test_approx_eq r[1] 1.0
println("solve x^3+log(2.0-x)-1=0 initial=50.0 in[0.0,100] with  attempts=",r[2])

@test_throws ErrorException Solver.callfzero(x->x^2+log(2.0-x)-1,-50.0,[-100.0,-10.0],100)

@test_throws DomainError Solver.callfzero(x->x^2+log(2.0-x)-1,500.0,[100.0,10000.0],100)

@time r=Solver.callfzero(y->(y)^3+log(y-1),50.0,[0.0,100.0],1000)
println("(y)^3+log(y-1) ->y=",r)

println("**** test opt ****")
f(x,y)=begin
        global output
        output && println("in f $x , $y");
        return x^3+log(y)
       end
g(x,y)=begin
        global output
        output && println("in g $x , $y");
        return log(x-y)
       end
indG=Array((Vector{Int}),0)
push!(indG,[1,2],[1,2])
# cant solve because it sends y>x and log(x-y) fails 
@test_throws DomainError Solver.callmin([f,g],indG,[0.0,0.0],[100.0,100.0],[10.0,1.0])

println("**** test fzero(f(x,fzero(g(x,y)))) instead of opt ****")
@time r=Solver.callfzero(x->f(x,Solver.callfzero(y->g(x,y),50.0,[0.0,100.0],1000)[1]),50.0,[0.0,100.0],1000)
println("use 2 fzeroz ->",r)

gg(vec_x)=g(vec_x...)
ff(vec_x)=f(vec_x...)
#output=true
println("**** test ffzero(f,g) instead of opt ****")
@time r=Solver.callffzero(ff,50.0,[0.0,100.0],gg,50.0,[0.0,100.0])
println("use ffzero ->",r)

println("**** test ffzero(g,f) instead of opt ****")
@time r=Solver.callffzero(gg,5.0,[0.0,10.0],ff,5.0,[0.0,10.0])
println("use ffzero ->",r)
