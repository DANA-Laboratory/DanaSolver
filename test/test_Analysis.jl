@test Solver.isconstantfactor([1.0,0.0,0.0,0.0,0.0,0.0])==true
function analysissystemofequations()
  println("**** analysis system of equations ****")
  exprs=Solver.Equations([
                           Solver.Equation(:(x+2*y+z^2)),
                           Solver.Equation(:(2*x-y*z)),
                           Solver.Equation(:(2*(z-2x)+3*log(y))),
                           Solver.Equation(:(x*y*z)),
                           Solver.Equation(:(x^y*log(z))),
                           Solver.Equation(:(x+y)),
                           Solver.Equation(:(3*z-y+x)),
                           Solver.Equation(:(5.0+6.0))
                           ])
  Solver.analysis(exprs)
  println(exprs.exs[6].ex,"-->",exprs.facts[1,:],'*',exprs.terms,"=0")
  println(exprs.exs[7].ex,"-->",exprs.facts[2,:],'*',exprs.terms,"=0")
  @test exprs.indexnonliexs==[1,2,3,4,5] #index of noliexprs
  #list of args
  for i in [1:8]
    @test exprs.exs[i].termall == [Set(["x","z","y"]),Set(["x","z","y"]),Set(["x","z","y"]),Set(["x","z","y"]),Set(["x","z","y"]),Set(["x","y"]),Set(["x","z","y"]),Set([])][i]
  #nonlinear terms of each nonlinear expr
    @test exprs.exs[i].termnonli == [Set(["z"]),Set(["z","y"]),Set(["y"]),Set(["x","z","y"]),Set(["x","z","y"]),Set{}(),Set{}(),Set{}()][i]
  end
  println("linear factors of nonlinear exprs:")
  for i in exprs.indexnonliexs
    ex=exprs.exs[i].ex
    facs=exprs.exs[i].factliinnonli
    println(ex,"-->",facs[2:end],'*',exprs.terms,"=0")
  end
end
function allsyms()
  println("**** get all syms in expr ****")
  syms=Set{String}()
  Solver.allsyms!(:(),syms)
  @test syms == Set{AbstractString}()
  syms=Set{String}()
  Solver.allsyms!(:(a+b),syms)
  @test syms == Set(AbstractString["b","a"]) 
  syms=Set{String}()
  Solver.allsyms!(:(-b*2*a+Sin(g+z^s-Cos(b+Log(3*-d)))),syms)
  @test syms == Set(AbstractString["b","a","d","g","s","z"]) 
end
