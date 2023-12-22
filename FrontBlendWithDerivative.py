import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import SR1

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
    foil.RampLower(0.1,0.1)

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
    i = res3000.LocalLowerBucketIndex()
    CaLow = res3000.Cl[i]
    errorLowerBucker = abs(CaLow - 0.3)
    weight = 40
    Error = Error + errorLowerBucker * weight
    print("Lower Bucket " + str(CaLow))

    ##Now at 2000
    foil = OptiWrapperGenerate(x,2000,alphaMin=4,alphaMax=10,open=True)
    res2000 = foil.ReadResults()
    #Upper Bucker
    i = res2000.LocalUpperBucketIndex()
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
        errorClMax = (1.5 - ClMax) * 40
    else:
        errorClMax = (1.5 - ClMax) * 5
    weight = 1
    Error = Error + errorClMax * weight
    print("ClMax " + str(ClMax))
    
    print("Error: " + str(Error))

    return Error

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


##This is very mehhh time for a custom optimizer
#xs = np.array([[6,0.1,8,0,0.3,0.7,3,0.15,2,0,0.3,0.7]])

#x0 = [6,0.1,8,0,0.3,0.7,3,0.15,2,0,0.3,0.7]

#bounds = Bounds([3, 0, 3, -5 ,0.3 ,0,0,0,-2,-5,0.3,0], [20, 1, 20, 5 , 1, 1.2, 6, 1, 5, 5, 1, 1])

#linear_constraint = LinearConstraint([[0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0]], [0],[np.inf])

#deleted min denominator
#res = minimize(Evaluate, x0, method='trust-constr',  jac=derivFoilEval, hess=SR1(min_denominator=0.05),constraints=[linear_constraint],options={'verbose': 0,"maxiter":25,"initial_tr_radius":1}, bounds=bounds)

#print(res)

#OptiWrapperGenerate(res.x,3000,open=True)

def singleChangeOptimizer(Eval,upperBounds,lowerBounds,iniGuess,iniStepSize,maxIter = 30):
    #Eval is function that we minimize
    #upper/lowerBounds are the bounds of the variables
    #iniGuess is the first guess for the Variables
    #iniStepSize is the initial step size
    
    n = len(upperBounds)
    iter = 0

    if(n != len(lowerBounds) or n != len(iniGuess) or n != len(iniStepSize)):
        print("Unput sizes dont match")
        return -1
    
    delta = np.ones(n)*99999 #delta is how much the error changed with the last change of this variable (normalized) with iniStepSize
    multiplier = np.ones(n)         #multiplier for the delta
    direction = np.ones(n)           #direction of the last step where error decreased

    stepSize = iniStepSize
    guess = iniGuess
    currentError = Eval(guess)

    targetStep = 1          #wir versuchen error um 1 zu ändern

    while(iter < maxIter):
        iter = iter+1

        print(currentError)

        #if(iter == n + 1):  #reselt multiplieres after first initial runn
         #   multiplier = np.ones(n)
        #for now we wana just do one after the other
        if(iter % n == 0):
            delta = np.ones(n)*99999 #delta is how much the error changed with the last change of this variable (normalized) with iniStepSize
            multiplier = np.ones(n)

        change = np.multiply(delta,multiplier)
        i = np.argmax(change,axis=None) #finde which change will probably be most effective
        #increse multiplier of all the others
        multiplier = np.add(multiplier,np.ones(n))
        multiplier[i] = 1

        #first step according to direction
        _guess = guess
        _guess[i] = _guess[i] + direction[i] * stepSize[i]

        #check constraints
        if(_guess[i] < lowerBounds[i]):
            _guess[i] = lowerBounds[i]
        elif(_guess[i] > upperBounds[i]):
            _guess[i] = upperBounds[i]

        _Error1 = Eval(_guess)
        if(_Error1 < currentError):
            delta[i] = abs(currentError - _Error1)
            currentError = _Error1
            guess = _guess
            #adapt stepsize so we move 1
            stepSize[i] = targetStep * stepSize[i] / delta[i]
            continue
        else:   #try other direction
            _guess[i] = _guess[i] - 2 * direction[i] * stepSize[i]
            _Error = Eval(_guess)
            if(_Error < currentError):
                delta[i] = abs(currentError - _Error)
                currentError = _Error
                direction[i] = -1 * direction[i]
                stepSize[i] = targetStep * stepSize[i] / delta[i]
                guess = _guess
                continue
            else:
                _delta = abs(_Error1 - _Error)

                stepSize[i] = targetStep * 2 * stepSize[i] / delta[i]

                delta[i] = 0   #punishment, since we dont improve
    print("Final best Error:" + str(currentError))
    return guess

#aLamT,startBlendT,LeT,dLeT,recTop,strengthRecTop,aLamB,startBlendB,LeB,dLeB,recBot,strengthRecBot
customGuess = [8/0.11,0.1,12,0,0.2,0.7,3,0.15,2,0,0.3,0.7]
upperBounds = [20, 0.35, 20, 5 , 0.9, 1.2, 6, 0.35, 4.5, 5, 1, 1]
lowerBounds = [3, 0.01, 5, -5 ,0.4 ,0,0,0.02,-2,-5,0.4,0]
iniGuess = [6,0.1,8,0,0.7,0.7,3,0.15,2,0,0.7,0.7]
iniStepSize = [1,0.1,1,0.5,0.2,0.2,1,0.1,0.5,0.5,0.2,0.2]
#bounds = Bounds([3, 0, 3, -5 ,0.3 ,0,0,0,-2,-5,0.3,0], [20, 1, 20, 5 , 1, 1.2, 6, 1, 5, 5, 1, 1])

#res = singleChangeOptimizer(Evaluate,upperBounds,lowerBounds,iniGuess,iniStepSize,maxIter=100)
#OptiWrapperGenerate(res,3000,open=True,num=7)
foil = OptiWrapperGenerate(iniGuess,3000,alphaMin=0,alphaMax=5,open=True)
self = foil.ReadResults()