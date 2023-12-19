import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil
import matplotlib.pyplot as plt

#first easy test
def FoilEval(aTop,aBottom,lengthTop,lengthBottom,MPRTop,MPRBottom) -> float:
    #create Foil
    foil = EpplerFoil.AirFoil("entwurf.dat","TF",1,_N = 60, _ZeroN = 31)
    foil1.BeginnFile()

    #Topside
    foil.WriteUpperC_Alpha([0],[aTop])
    foil.MPRUpper(lengthTop,0.98,-1,MPRTop,mode=2)
    foil.RampUpper(0.1,0.1)

    #Bottomside
    foil1.WriteLowerC_Alpha([1],[aBottom])
    foil1.MPRLower(lengthBottom,0.98,0.23,MPRBottom)
    foil1.RampLower(0.1,0.1)

    #Finalize
    foil.writeMPR("upper",0.4,0)