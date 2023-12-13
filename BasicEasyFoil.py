import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil

##create a laminar airfoil from basic variables

##Bottom side
lamBottom = 0.3 #Ca at which the bottom surface should be laminar
lamLengthBottom = 0.7 #position where MPR starts at bottom -> lower end of bucket
#Ramp Bottom Side against lam seperation
blendingLamBottom = 0.5 #reach of ramp into lam area -> increase when laminar separation occurs at MPR
blendingMPRBottom = 0.1 #ramp into MPR -> to combat laminar separation 
#MPR bottom side 
shapeMPRBottom = 0.23
strengthMPRBottom = 0.7

##Top Side
lamTop = 0.8 #Target Ca for laminar flow at Top -> upper end of bucket
lamLengthTop = 0.7 #laminar length
#Ramp Top Side against lam seperation
blendingLamTop = 0.5 #reach of ramp into lam area -> increase when laminar separation occurs at MPR
blendingMPRTop = 0.1 #ramp into MPR -> to combat laminar separation 
#MPR Top side 
shapeMPRTop = -1
strengthMPRTop = 0.7
#alpah* increase against suction peak -> increases alpha* smothly (according to a polynom) towards the LE to combat suction peaks
alphaLE = 0 #the alpha* we will have at the LE
startIncrease = 0.3  #the c at which the increase starts
numberPoints = 8    #amount of points to use

adaptSide = "upper" #which side to iterate for MPR choices are "upper", "lower", "mixed", "none"


##GenerateFile
foil1 = EpplerFoil.AirFoil("entwurf.dat","TF",420)
foil1.BeginnFile()

#generate upper side alpha*
#fit function to combat LE
#func is a*x**3 + b*x**2 + c*x + d

opti = asb.Opti() #erstellt optimizer environment

a = opti.variable(init_guess=1)
b = opti.variable(init_guess=1)
c = opti.variable(init_guess=1)
d = opti.variable(init_guess=1)

#function to interpolate:
def f(x,a,b,c,d):
    return a*x**3 + b*x**2 + c*x + d
def df(x,a,b,c,d):
    return 3*a*x**2 + 2*b*x + c
def ddf(x,a,b,c,d):
    return 6*a*x + 2*b

opti.subject_to(f(startIncrease,a,b,c,d) == (lamTop/0.11))
opti.subject_to(f(0,a,b,c,d) == alphaLE)
opti.subject_to(df(startIncrease,a,b,c,d) == 0)
opti.subject_to(ddf(startIncrease,a,b,c,d) == 0)

sol = opti.solve()
a = sol.value(a)
b = sol.value(b)
c = sol.value(c)
d = sol.value(d)

xarr = np.linspace(startIncrease,0,num=numberPoints)
yarr = np.zeros(numberPoints)
for i in range(0,numberPoints):
    x = xarr[i]
    y = a*x**3 + b*x**2 + c*x + d
    yarr[i] = y

foil1.WriteUpperC_Alpha(xarr,yarr)
foil1.MPRUpper(lamLengthTop,0.98,shapeMPRTop,strengthMPRTop)
foil1.RampUpper(blendingLamTop,blendingMPRTop)

#Bottom Side
foil1.WriteLowerC_Alpha([1],[(lamBottom/0.11)])
foil1.MPRLower(lamLengthBottom,0.98,shapeMPRBottom,strengthMPRBottom)
foil1.RampLower(blendingLamBottom,blendingMPRBottom)
foil1.writeMPR(adaptSide,0.4,0)


#Analysis
foil1.inviscidCalc_byIncrements(0,3,5)
foil1.naturalTransitionCalc(3000)
foil1.visousCalc(0,15,15)
foil1.CloseFile()
foil1.ExecuteEppler()





