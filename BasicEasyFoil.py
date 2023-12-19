import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil
import matplotlib.pyplot as plt

##create a laminar airfoil from basic variables

##Bottom side
lamBottom = 0.4 #Ca at which the bottom surface should be laminar
lamLengthBottom = 0.70 #position where MPR starts at bottom -> lower end of bucket
#Ramp Bottom Side against lam seperation
blendingLamBottom = 0.1 #reach of ramp into lam area -> increase when laminar separation occurs at MPR
blendingMPRBottom = 0.1 #ramp into MPR -> to combat laminar separation 
#MPR bottom side 
shapeMPRBottom = 0.23
strengthMPRBottom = 0.7
#anti suction peak bottom side
cStart = 0.01
cEnde = 0.1
alphaBE = lamBottom/0.11 #here a smaller number combats suction peaks
numBottom = 5 #amount of points

##Top Side
lamTop = 0.8 #Target Ca for laminar flow at Top -> upper end of bucket
lamLengthTop = 0.6 #laminar length
#Ramp Top Side against lam seperation
blendingLamTop = 0.1 #reach of ramp into lam area -> increase when laminar separation occurs at MPR
blendingMPRTop = 0.1 #ramp into MPR -> to combat laminar separation 
#MPR Top side 
shapeMPRTop = -1
strengthMPRTop = 0.7
#alpah* increase against suction peak -> increases alpha* smothly (according to a polynom) towards the LE to combat suction peaks
alphaLE = lamTop/0.11 #the alpha* we will have at the LE
startIncrease = 0.2  #the c at which the increase starts
numberPoints = 4    #amount of points to use

adaptSide = "upper" #which side to iterate for MPR choices are "upper", "lower", "mixed", "none"


##GenerateFile
foil1 = EpplerFoil.AirFoil("entwurf.dat","TF",427)
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

#Bottom Anti suction peak
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

opti.subject_to(f(cEnde,a,b,c,d) == (lamBottom/0.11))
opti.subject_to(f(cStart,a,b,c,d) == alphaBE)
opti.subject_to(df(cEnde,a,b,c,d) == 0)
opti.subject_to(ddf(cEnde,a,b,c,d) == 0)

sol = opti.solve()
a = sol.value(a)
b = sol.value(b)
c = sol.value(c)
d = sol.value(d)

xarr = np.linspace(cStart,cEnde,num=numBottom)
_x = np.zeros(numBottom+1)
yarr = np.zeros(numBottom+1)
for i in range(0,numBottom):
    x = xarr[i]
    _x[i] = x
    y = a*x**3 + b*x**2 + c*x + d
    yarr[i] = y
_x[-1] = 1
yarr[-1] = (lamBottom/0.11)

#Bottom Side
foil1.WriteLowerC_Alpha(_x,yarr)
foil1.MPRLower(lamLengthBottom,0.98,shapeMPRBottom,strengthMPRBottom)
foil1.RampLower(blendingLamBottom,blendingMPRBottom)
foil1.writeMPR(adaptSide,0.4,0)


#Analysis
foil1.inviscidCalc_byIncrements(0,3,7)
foil1.forcedTransitionCalc(3000,100,75)
#foil1.forcedTransitionCalc(2000,100,75)
#foil1.forcedTransitionCalc(1000,100,75)
foil1.visousCalc(0,12,10)
foil1.CloseFile()
foil1.ExecuteAndOpen()
#foil1.Execute()
#

calcRes = foil1.ReadResults()
iUp = calcRes.UpperBucketIndex()
iLow = calcRes.LowerBucketIndex()
print(calcRes.Cd)