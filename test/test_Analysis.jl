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
function testanalysis()
  equation=[:(teta - acos(r / q ^ 1.5)),:(Z1 - (-2 * sqrt(q) * cos(teta / 3) - beta / 3)),:(Z2 - (-2 * sqrt(q) * cos((6.283185307179586 + teta) / 3) - beta / 3)),:(Z3 - (-2 * sqrt(q) * cos((12.566370614359172 + teta) / 3) - beta / 3)),:(1.0 - ((Tr / 0.7499770413536668) / 0.31115542517055556 - (4.722638593459387alpha) / 1.4375344374424004)),:(Z - 0.31115542517055556 / Tr),:(A - (0.0968176986130692 * (4.722638593459387alpha)) / Tr ^ 2),:(alpha - (1 + 0.8001649902192 * (1 - sqrt(Tr))) ^ 2),:(B - 0.077796 / Tr),:(beta - (B - 1)),:(gama - ((A - 3 * B ^ 2) - 2B)),:(delta - ((B ^ 3 + B ^ 2) - A * B)),:(q - (beta * beta - 3gama) / 9),:(r - ((2 * beta ^ 3 - 9 * beta * gama) + 27delta) / 54),:(s_Dep - (0.9670056591193716 * (1.8001649902192 / sqrt(Tr) - 0.8001649902192) - log(Z - B))),:(h_Dep - ((1 - Z) + 2.1755134930530544 * sqrt(alpha))),:(g_Dep - ((1 - Z) + (0.8546225611932206alpha) / (0.7071717180445836Tr) + log(0.7499770413536668Z)))]

  println(Solver.analysis([:(x^2)]))
  println(Solver.analysis([:(1/x^2)]))
  println(Solver.analysis([:(x/2)]))
  println(Solver.analysis([:(1*x^2)]))
  return
  for eq in equation
    println(eq)
    re=Solver.analysis([eq])
    println(re[end])
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
