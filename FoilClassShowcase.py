#test import
import EpplerFoil
import numpy as np

##Setup

#first create a foil handle, set the file it will be saved in "entwurf.dat"
#and name the foil this one will be called "TF 187" and have 60 Points
#additionally since the foil isnt a true circle we add a guess of what N corresponds to the LE
#-> if this guess it to high the Eppler Code will crash 
#_N and _ZeroN can be omitted and are set to standard values 60 and 31
foil1 = EpplerFoil.AirFoil("entwurf.dat","TF",187,_N = 60, _ZeroN = 31)

#Now a file is created and the REMO lines are written, allways call this before the remaining commands
foil1.BeginnFile()  #start writing to file


##Defining Upper Surface

#First we set our alpha*s, here the coordinates are relative to chord length and NOT in polar coordinates
#here we want from the TE to c=0.3 a alpha* of 8.3 and from 0.3 to 0.1 a alpha* of 10 -> we work from the back to the front
#0 can only be used for the LE (its weird but it is what it is) and must be pressent, so if its not part of the input
#it is automatically added. A warning will be displayed if 0 point is added to remind you to only use it for the LE
foil1.WriteUpperC_Alpha([0.3,0.1],[8.3,10])

#Now we create the MPR, first the chord lenght where it should start (here c = 0.6),
#then we say at what length the closing curve should be added (here c = 0.98 <- is a good standard value).
#afterwards we set the shape of the curve (-1 see eppler) and strength(0.7 again see eppler)
#the mode can be specified, but defaults to 2 if not -> see eppler user guide
foil1.MPRUpper(0.65,0.98,-1,0.7,mode=2)

#Now we specify the shape of the ramp to smoothly go into the MPR,
#coordinates here are a deltaC -> we want the ramp to start 0.2c infront of the MPR and end 0.05c after it
foil1.RampUpper(0.2,0.05)

##Defining Lower Surface

#follows the same rules however the c's must now be increasing -> from front to back
foil1.WriteLowerC_Alpha([0.98],[2.7])
foil1.MPRLower(0.7,0.98,0.23,0.7)
foil1.RampLower(0.2,0.05)

##Create Airfoil
#This creates the Airfoil and finalizes the MPR line, first we select the mode
#"upper" -> adapt upper alphas, "lower" -> adapt lower side alphas, "mixed" -> adapt both
#all the other modes that are in the user guide can just be input as number
#Next up we have the sum of the Ks, 0.4 seems like a good value, mainly influences the shape of the closing curve
#and finaly the tolerance is specified, defaults to 0
foil1.writeMPR("upper",0.4,0)

##Calculations
#here are multiple ways to perform a velocity calculation
#foil1.inviscidCalc_Custom([0,3,6,8]) -> enter all the aoa's (respective to 0 lift) yourself as an array
#foil1.invisidCalc(0,15,5) -> aoa's from 0 to 15 degree in 5 steps
foil1.inviscidCalc_byIncrements(0,3,5) #-> start at an aoa of 0 and increase by 3 each time for a total of 5 points

#Now to execute the BL calculations, here the calculations are added for each RE and transition you want and then its executed
#adds a natural Transition Calculation for RE = 3.000.000
foil1.naturalTransitionCalc(3000)
#adds a calculation with forced transition for RE = 2.000.000 with Uper transition at 100% chord(aka natural transition) and 75% chord at the bottom side
foil1.forcedTransitionCalc(2000,100,75)
#now we execture the calculation from an aoa of 0 degree to 15 degree in 15 steps
foil1.visousCalc(0,15,15)
#Cause I like these settings for the BL calculation the BL values are written to the .l file -> you can check where seperation, transition, reattachment etc. happens
#And additionally triangles are added to the plot that symbolize if seperation occurs at the upper or lower surface

#Now we add the END line and save the file now you can use it in the eppler code
foil1.CloseFile()

#Now run the Eppler Code with the created file
foil1.ExecuteEppler()

