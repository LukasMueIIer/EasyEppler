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
    foil = EpplerFoil.AirFoil("entwurf.dat","AO",10,_N = 60, _ZeroN = 31)
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
    foil.visousCalc(0,17,17)
    foil.CloseFile()
    if(open):
        foil.ExecuteAndOpen()
    else:
        foil.Execute()
    return foil

def OptiWrapperGenerate(x,RE,open=False,num=7):
    return GenerateFoil(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],RE=RE,open=open,num=num)

OptiWrapperGenerate([6,0.1,8,0,0.3,0.7,3,0.15,2,0,0.3,0.7],RE=3000,open=True,num=7)


def Evaluate(x):
    ##We evaluate for multiple RE numbers starting with 3000
    foil = OptiWrapperGenerate(x,3000)
    res3000 = foil.ReadResults()
    #Thickness
    targThick = 15
    errorThick = abs(targThick - res3000.thickness)
    weight = 1
    Error = errorThick * weight
    print("Thickness:" + str(res3000.thickness))
    #Bucket at ca = 0.3
    i = res3000.LowerBucketIndex()
    CaLow = res3000.Cl[i]
    errorLowerBucker = abs(CaLow - 0.3)
    weight = 20
    Error = Error + errorLowerBucker * weight
    print("Lower Bucket " + str(CaLow))

    ##Now at 2000
    foil = OptiWrapperGenerate(x,2000)
    res2000 = foil.ReadResults()
    #Upper Bucker
    i = res2000.UpperBucketIndex()
    CaHigh = res2000.Cl[i]
    errorUpperBucker = abs(CaHigh - 0.8)
    weight = 20
    Error = Error + errorUpperBucker * weight
    print("Upper Bucket " + str(CaHigh))

    #Now at 1000
    foil = OptiWrapperGenerate(x,1000)
    res1000 = foil.ReadResults()
    #high lift
    ClMax = res1000.MaxCl()
    if(ClMax < 1.5): #Non Linear punishment for underperforming
        errorClMax = (1.5 - ClMax) * 6
    else:
        errorClMax = (1.5 - ClMax) 
    weight = 10
    Error = Error + errorClMax * weight
    print("ClMax " + str(ClMax))
    
    print("Error: " + str(Error))

    return Error