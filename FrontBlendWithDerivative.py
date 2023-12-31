import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import SR1
import casadi as cas

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
 

    sol = opti.solve(verbose=False)
    a = sol.value(a)
    b = sol.value(b)
    c = sol.value(c)
    d = sol.value(d)
    e = sol.value(e)

    res = np.zeros(len(x))
    for i in range(0,len(x)):
        res[i] = f(x[i],a,b,c,d,e)

    return res

def GenerateFoil(aLamT,startBlendT,LeT,dLeT,recTop,strengthRecTop,aLamB,startBlendB,LeB,dLeB,recBot,strengthRecBot,RE=3000,open=False,num=7,alphaMin=0,alphaMax=17) -> EpplerFoil.AirFoil:
    foil = EpplerFoil.AirFoil("entwurf.dat","AO",12,_N = 60, _ZeroN = 31)
    foil.BeginnFile()

    #Upper Side
    x = np.linspace(startBlendT,0,num=num)
    y = funcFit(startBlendT,aLamT,0,0,0,LeT,dLeT,x)
    foil.WriteUpperC_Alpha(x,y)
    foil.MPRUpper(recTop,0.98,-1,strengthRecTop,mode=2)
    foil.RampUpper(0.1,0.1)                         #TODO ADD RAMP

    #LowerSide
    x = np.linspace(0.01,startBlendB,num=num)
    y = funcFit(startBlendB,aLamB,0,0,0,LeB,dLeB,x)
    foil.WriteLowerC_Alpha(x,y)
    foil.MPRLower(recBot,0.98,0.23,strengthRecBot,mode=2)
    foil.RampLower(0.05,0.05)

    #Finalize
    foil.writeMPR("upper",0.4,0)

    #Calculations
    foil.inviscidCalc_byIncrements(0,3,5)
    foil.naturalTransitionCalc(RE)
    foil.visousCalc(alphaMin,alphaMax,15)
    foil.CloseFile()
    if(open):
        foil.ExecuteAndOpen()
    else:
        foil.Execute()
    return foil

def OptiWrapperGenerate(x,RE,open=False,num=7,alphaMin=0,alphaMax=17):
    print(x)
    return GenerateFoil(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],RE=RE,open=open,num=num,alphaMin=alphaMin,alphaMax=alphaMax)

def Evaluate(x):
    ##We evaluate for multiple RE numbers starting with 3000
    foil = OptiWrapperGenerate(x,3000,alphaMin=0,alphaMax=5,open=True)
    res3000 = foil.ReadResults()
    #Thickness
    targThick = 15
    errorThick = abs(targThick - res3000.thickness)
    weight = 10
    Error = errorThick * weight
    print("Thickness:" + str(res3000.thickness))
    #Bucket at ca = 0.3
    i = res3000.LowerBucketIndex()
    CaLow = res3000.Cl[i]
    errorLowerBucker = abs(CaLow - 0.3)
    weight = 40
    Error = Error + errorLowerBucker * weight
    print("Lower Bucket " + str(CaLow))

    ##Now at 2000
    foil = OptiWrapperGenerate(x,2000,alphaMin=4,alphaMax=10,open=True)
    res2000 = foil.ReadResults()
    #Upper Bucker
    i = res2000.UpperBucketIndex()
    CaHigh = res2000.Cl[i]
    errorUpperBucker = abs(CaHigh - 0.8)
    weight = 40
    Error = Error + errorUpperBucker * weight
    print("Upper Bucket " + str(CaHigh))

    #Now at 1000
    foil = OptiWrapperGenerate(x,1000,alphaMin=6,alphaMax=17)
    res1000 = foil.ReadResults()
    #high lift
    ClMax = res1000.MaxCl()
    if(ClMax < 1.5): #Non Linear punishment for underperforming
        errorClMax = (1.5 - ClMax) * 10
    else:
        errorClMax = (1.5 - ClMax) * 5
    weight = 1
    Error = Error + errorClMax * weight
    print("ClMax " + str(ClMax))
    
    print("Error: " + str(Error))

    return Error

def EvaluateVec(x): #returns the current evaluation vector
    print("EVALUATE")
    print(x)
    res = np.zeros(4) #our vector has 4 entries
    ##We evaluate for multiple RE numbers starting with 3000
    foil = OptiWrapperGenerate(x,3000,alphaMin=0,alphaMax=5.5,open=False)
    res3000 = foil.ReadResults()
    #Thickness
    print("Thickness:" + str(res3000.thickness))
    res[0] = res3000.thickness
    #Bucket at ca = 0.3
    i = res3000.LowerBucketIndex()
    CaLow = res3000.Cl[i]
    res[1] = CaLow
    print("Lower Bucket " + str(CaLow))

    ##Now at 2000
    foil = OptiWrapperGenerate(x,2000,alphaMin=4,alphaMax=10,open=False)
    res2000 = foil.ReadResults()
    #Upper Bucker
    i = res2000.UpperBucketIndex()
    CaHigh = res2000.Cl[i]
    res[2] = CaHigh
    print("Upper Bucket " + str(CaHigh))

    #Now at 1000
    foil = OptiWrapperGenerate(x,1000,alphaMin=6,alphaMax=17)
    res1000 = foil.ReadResults()
    #high lift
    ClMax = res1000.MaxCl()
    res[3] = ClMax 
    print("ClMax " + str(ClMax))
    
    return res

def derivFoilEval(x,delta=[1,0.03,0.5,0.2,0.05,0.1,1,0.03,0.05,0.2,0.05,0.1]) -> float: #scyPy takes way to small steps that get removed in eppler rounding
    res = np.zeros(len(x))
    for i in range(0,len(x)): #we use centrall diff approach
        _x = x
        _x[i] = _x[i] + 0.5 * delta[i]
        UpScore = Evaluate(_x)
        _x[i] = _x[i] - 0.5 * delta[i]
        LowScore = Evaluate(_x)
        res[i] = (UpScore - LowScore)/delta[i]
    print(res)
    return res

def VectorwiseOptimizer(Eval,target,lowerBounds,upperBounds,iniGuess,relaxation,gradientStepSize,maxIterations = 10): #BOUNDS MUST FULLFILL REQUIREMENTS FROM EPPLER CODE
    currEval = Eval(iniGuess)
    currGuess = iniGuess
    iter = True
    iterCount = 0
    n = len(iniGuess)
    while(iter):
        iter = iterCount + 1
        if(iter >= maxIterations):
            iter = False
            continue
        #Build gradient
        gradMatrix = np.zeros((len(target),n)).T
        for i in range(0,n):    #for performance now one sided derivative
            derivGuess = np.copy(currGuess)
            derivGuess[i] = derivGuess[i] + gradientStepSize[i]

            if(derivGuess[i] < lowerBounds[i]): #check boundaries
                derivGuess[i] = lowerBounds[i]
            elif(derivGuess[i] > upperBounds[i]):
                derivGuess[i] = upperBounds[i]

            #calculate derivative
            deltaX = derivGuess[i] - currGuess[i]
            print("GUESSSS")
            print(derivGuess)
            deltaY = np.subtract(Eval(derivGuess),currEval)
            _deriv = (1/deltaX) * deltaY

            #write into Matrix
            gradMatrix[i] = _deriv

        #now solve with optimizer
        opt = asb.Opti()
        x = opt.variable(init_guess=np.zeros(n)) #our vector
        opt.subject_to((currGuess + x) > lowerBounds)
        opt.subject_to((currGuess + x) < upperBounds)
        #opt.subject_to(np.abs(x) < gradientStepSize)    #This sucks underrelaxation is the way to go
        opt.minimize(np.sum((currEval + cas.mtimes(x.T,gradMatrix).T - target)**2)) #we minimize insted of subject
        res = opt.solve()
        delta = res.value(x)

        #now do the step
        currGuess = np.add(currGuess,relaxation*delta)
        currEval = Eval(currGuess)
        print(currEval)

    print("Final Error")
    print(currEval)
    return currGuess

x0 = [6,0.1,8,0,0.7,0.7,3,0.15,2,0,0.75,0.7]
#aLamT,startBlendT,LeT,dLeT,recTop,strengthRecTop,aLamB,startBlendB,LeB,dLeB,recBot,strengthRecBot
lowerBounds = [3, 0.02, 3, -5 ,0.35 ,0,0,0.02,-2,-5,0.74,0]
upperBounds = [10, 0.3, 20, 5 , 0.9, 1, 6, 0.3, 5, 5, 0.76, 1]
delta = [0.8,0.03,0.8,0.05,0.1,0.1,0.8,0.03,0.8,0.05,0.1,0.1]
res = VectorwiseOptimizer(EvaluateVec,[15,0.3,0.8,1.5],lowerBounds,upperBounds,x0,0.05,delta,maxIterations=200)

#OptiWrapperGenerate(x0,3000,open=True)
print("TEST")
