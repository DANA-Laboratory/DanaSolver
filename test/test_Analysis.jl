@test Solver.isconstantfactor([1.0,0.0,0.0,0.0,0.0,0.0])==true
function analysissystemofequations()
  println("**** analysis system of equations ****")
  exprs=[:(x+2*y+z^2),:(2*x-y*z),:(2*(z-2x)+3*log(y)),:(x*y*z),:(x^y*log(z)),:(x+y),:(3*z-y+x),:(5.0+6.0)]
  re=Solver.analysis(exprs)
  println(exprs[6],"-->",re[1][1,:],'*',re[2],"=0")
  println(exprs[7],"-->",re[1][2,:],'*',re[2],"=0")
  @test length(re[4])==8 #list of each expr args
  @test re[3]==[1,2,3,4,5] #index of noliexprs
  #list of args
  @test re[4] == [Set(["x","z","y"]),Set(["x","z","y"]),Set(["x","z","y"]),Set(["x","z","y"]),Set(["x","z","y"]),Set(["x","y"]),Set(["x","z","y"]),Set([])]
  #nonlinear terms of each nonlinear expr
  @test re[5] == [Set(["z"]),Set(["z","y"]),Set(["y"]),Set(["x","z","y"]),Set(["x","z","y"]),Set{}(),Set{}(),Set{}()]
  println ("linear factors of nonlinear exprs:")
  nonliexprs=exprs[re[3]]
  nonliexprsfacs=(re[6])[re[3]]
  for i in [1:length(nonliexprs)]
    ex=nonliexprs[i]
    facs=nonliexprsfacs[i]
    println(ex,"-->",facs[2:end],'*',re[2],"=0")
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
