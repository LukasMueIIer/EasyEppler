import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import SR1


#first easy test
def FoilEval(aTop,aBottom,lengthTop,lengthBottom,MPRTop,MPRBottom) -> float:
    #create Foil
    foil = EpplerFoil.AirFoil("entwurf.dat","TF",1,_N = 60, _ZeroN = 31)
    foil.BeginnFile()

    #Topside
    foil.WriteUpperC_Alpha([0.1,0],[aTop,15])
    foil.MPRUpper(lengthTop,0.98,-1,MPRTop,mode=2)
    foil.RampUpper(0.1,0.1)

    #Bottomside
    foil.WriteLowerC_Alpha([1],[aBottom])
    foil.MPRLower(lengthBottom,0.98,0.23,MPRBottom)
    foil.RampLower(0.1,0.1)

    #Finalize
    foil.writeMPR("upper",0.4,0)

    #Calculations
    foil.inviscidCalc_byIncrements(0,3,5)
    foil.naturalTransitionCalc(3000)
    foil.visousCalc(0,15,15)
    foil.CloseFile()
    foil.Execute()
    #foil.ExecuteAndOpen()
    calcRes = foil.ReadResults()

    #Evaluate
    targThick = 15
    dThick = abs(targThick - calcRes.thickness)
    Cd1 = calcRes.CdFromCl(0.3) 
    Cd2 = calcRes.CdFromCl(0.6)
    Clm = calcRes.MaxCl()
    score = dThick + 1000 * Cd1 + 1000 * Cd2 - Clm #beinhaltet Gewichtung für Optimierung
    print("Thickness:" + str(calcRes.thickness)+"% Cd1:" + str(Cd1) +" Cd2:" + str(Cd2) +" ClMax:" + str(Clm) + " Score:" + str(score))
    return score

def FoilShow(aTop,aBottom,lengthTop,lengthBottom,MPRTop,MPRBottom) -> float:
    #create Foil
    foil = EpplerFoil.AirFoil("entwurf.dat","TF",1,_N = 60, _ZeroN = 31)
    foil.BeginnFile()

    #Topside
    foil.WriteUpperC_Alpha([0.1,0],[aTop,15])
    foil.MPRUpper(lengthTop,0.98,-1,MPRTop,mode=2)
    foil.RampUpper(0.1,0.1)

    #Bottomside
    foil.WriteLowerC_Alpha([1],[aBottom])
    foil.MPRLower(lengthBottom,0.98,0.23,MPRBottom)
    foil.RampLower(0.1,0.1)

    #Finalize
    foil.writeMPR("upper",0.4,0)

    #Calculations
    foil.inviscidCalc_byIncrements(0,3,5)
    foil.naturalTransitionCalc(3000)
    foil.visousCalc(0,15,15)
    foil.CloseFile()
    #foil.Execute()
    foil.ExecuteAndOpen()
    calcRes = foil.ReadResults()

    #Evaluate
    targThick = 15
    dThick = abs(targThick - calcRes.thickness)
    Cd1 = calcRes.CdFromCl(0.3) 
    Cd2 = calcRes.CdFromCl(0.6)
    Clm = calcRes.MaxCl()
    score = dThick + 1000 * Cd1 + 1000 * Cd2 - Clm #beinhaltet Gewichtung für Optimierung
    print("Thickness:" + str(calcRes.thickness)+"% Cd1:" + str(Cd1) +" Cd2:" + str(Cd2) +" ClMax:" + str(Clm) + " Score:" + str(score))
    return score

def ScyFoilEval(x) -> float:
    return FoilEval(x[0],x[1],x[2],x[3],x[4],x[5])

def FoilFinal(x) -> float:
    return FoilShow(x[0],x[1],x[2],x[3],x[4],x[5])

def derivFoilEval(x,delta) -> float: #scyPy takes way to small steps that get removed in eppler rounding
    res = np.zeros(len(x))
    for i in range(0,len(x)): #we use centrall diff approach
        _x = x
        _x[i] = _x[i] + 0.5 * delta[i]
        UpScore = ScyFoilEval(_x)
        _x[i] = _x[i] - 0.5 * delta[i]
        LowScore = ScyFoilEval(_x)
        res[i] = (UpScore - LowScore)/delta[i]
    print(res)
    return res


     

#SciPy Optmizer

#(aTop,aBottom,lengthTop,lengthBottom,MPRTop,MPRBottom)
x0 = [5,1,0.7,0.7,0.7,0.7]


def use(x):
    delta = [1,1,0.1,0.1,0.1,0.1]
    return derivFoilEval(x,delta)

bounds = Bounds([0, 0, 0.15, 0.15 ,-1 ,-1], [20, 20, 0.9, 0.9 , 1, 1])
linear_constraint = LinearConstraint([[1, -1,0,0,0,0]], [0],[np.inf])

res = minimize(ScyFoilEval, x0, method='trust-constr',  jac=use, hess=SR1(min_denominator=0.05),constraints=[linear_constraint],options={'verbose': 0,"maxiter":25,"initial_tr_radius":0.1}, bounds=bounds)
print(res)
FoilFinal(res.x)