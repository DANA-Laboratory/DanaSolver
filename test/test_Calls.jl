
@test_approx_eq Solver.callfzero(x->x^3-1,0.0,[0.0,100.0])[1] 1.0

@time r=Solver.callfzero(x->x^2-1,200.0,[0.0,500.0])
println("attempts=$r")
@test_approx_eq Solver.callfzero(x->x^2-1,0.0,[0.0,100.0])[1] 1.0

@time r=Solver.callfzero(x->x^2+log(2.0-x)-1,50.0,[0.0,100.0],1000)
@test_approx_eq r[1] 1.0
println("attempts=$r")

@test_throws ErrorException Solver.callfzero(x->x^2+log(2.0-x)-1,-50.0,[-100.0,-10.0],100)

@test_throws DomainError Solver.callfzero(x->x^2+log(2.0-x)-1,500.0,[100.0,10000.0],100)

