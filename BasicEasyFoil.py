import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil

##create a laminar airfoil from basic variables

##Bottom side
lamBottom = 0.3 #Ca at which the bottom surface should be laminar
lamLengthBottom = 0.7 #position where MPR starts at bottom -> lower end of bucket
#Ramp Bottom Side against lam seperation
blendingLamBottom = 0 #reach of ramp into lam area -> increase when laminar separation occurs at MPR
blendingMPRBottom = 0 #ramp into MPR -> to combat laminar separation 
#MPR bottom side 
shapeMPRBottom = 0.23
strengthMPRBottom = 0.7

##Top Side
lamTop = 0.8 #Target Ca for laminar flow at Top -> upper end of bucket
lamLengthTop = 0.7 #laminar length
#Ramp Top Side against lam seperation
blendingLamTop = 0 #reach of ramp into lam area -> increase when laminar separation occurs at MPR
blendingMPRTop = 0 #ramp into MPR -> to combat laminar separation 
#MPR Top side 
shapeMPRTop = -1
strengthMPRTop = 0.7
#alpah* increase against suction peak -> increases alpha* smothly (according to a polynom) towards the LE to combat suction peaks
alphaLE = 0.8 / 0.11 #the alpha* we will have at the LE
startIncrease = 0.2  #the c at which the increase starts

##GenerateFile
foil = EpplerFoil.AirFoil("entwurf.dat","TF",420)
foil1.BeginnFile()
