import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil
import matplotlib.pyplot as plt

def funcFit(x0,Fx0,dFx0,ddFx0,x1,Fx1,dFx1,x):
    opti = asb.Opti() #erstellt optimizer environment

    a = opti.variable(init_guess=1)
    b = opti.variable(init_guess=1)
    c = opti.variable(init_guess=1)
    d = opti.variable(init_guess=1)
    e = opti.variable(init_guess=1)

    #function to interpolate:
    def f(x,a,b,c,d,e):
        return a*x**4 + b*x**3 + c*x**2 + d*x + e
    def df(x,a,b,c,d,e):
        return 4*a*x**3 + 3*b*x**2 + 2*c*x + d
    def ddf(x,a,b,c,d,e):
        return 12*a*x**2 + 6*b*x + 2*c

    opti.subject_to(f(x0,a,b,c,d,e) == (Fx0))
    opti.subject_to(df(x0,a,b,c,d,e) == (dFx0))
    opti.subject_to(ddf(x0,a,b,c,d,e) == (ddFx0))
    opti.subject_to(f(x1,a,b,c,d,e) == (Fx1))
    opti.subject_to(df(x1,a,b,c,d,e) == (dFx1))
 

    sol = opti.solve()
    a = sol.value(a)
    b = sol.value(b)
    c = sol.value(c)
    d = sol.value(d)
    e = sol.value(e)

    res = np.zeros(len(x))
    for i in range(0,len(x)):
        res[i] = f(x[i],a,b,c,d,e)

    return res

def GenerateFoil(aLamT,startBlendT,LeT,dLeT,recTop,strengthRecTop,aLamB,startBlendB,LeB,dLeB,recBot,strengthRecBot,RE=3000,open=False,num=7) -> EpplerFoil.AirFoil:
    foil = EpplerFoil.AirFoil("entwurf.dat","AO",5,_N = 60, _ZeroN = 31)
    foil.BeginnFile()

    #Upper Side
    x = np.linspace(startBlendT,0,num=num)
    y = funcFit(startBlendT,aLamT,0,0,0,LeT,dLeT,x)
    foil.WriteUpperC_Alpha(x,y)
    foil.MPRUpper(recTop,0.98,-1,strengthRecTop,mode=2)
    foil.RampUpper(0.1,0.1)                         #TODO ADD RAMP

    #LowerSide
    x = np.linspace(0.05,startBlendB,num=num)
    y = funcFit(startBlendB,aLamB,0,0,0,LeB,dLeB,x)
    foil.WriteLowerC_Alpha(x,y)
    foil.MPRLower(recBot,0.98,0.23,strengthRecBot,mode=2)
    foil.RampLower(0.1,0.1)

    #Finalize
    foil.writeMPR("upper",0.4,0)

    #Calculations
    foil.inviscidCalc_byIncrements(0,3,5)
    foil.naturalTransitionCalc(RE)
    foil.visousCalc(0,15,15)
    foil.CloseFile()
    if(open):
        foil.ExecuteAndOpen()
    else:
        foil.Execute()
    return foil

def OptiWrapperGenerate(x,RE,open=False,num=7):
    GenerateFoil(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],RE=RE,open=open,num=num)

OptiWrapperGenerate([6,0.1,8,0,0.3,0.7,3,0.15,2,0,0.3,0.7],RE=3000,open=True,num=7)