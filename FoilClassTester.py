#test import
import EpplerFoil
import numpy as np


foil1 = EpplerFoil.AirFoil("entwurf.dat","LM",187) #create Foil Handle
 
foil1.BeginnFile()  #start writing to file

#Set Upper Side
foil1.WriteUpperC_Alpha([0.3,0],[8.3,10])
foil1.MPRUpper(0.7,0.98,-1,0.7)
foil1.RampUpper(0.2,0.05)

#Set Lower Side
foil1.WriteLowerC_Alpha([0.98],[2.7])
foil1.MPRLower(0.7,0.98,0.23,0.7)
foil1.RampLower(0.2,0.05)

#Create MPR
foil1.writeMPR("upper",0.4,0)

#Calculations
#foil1.invisidCalc(0,15,5)
foil1.inviscidCalc_byIncrements(0,3,5)


foil1.CloseFile()


print("Test")
