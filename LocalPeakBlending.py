import aerosandbox as asb
import aerosandbox.numpy as np
import EpplerFoil

##create a laminar airfoil from basic variables
##here the anti suction peak function is modified
##so far not very good results tbh

##Bottom side
lamBottom = 0.3 #Ca at which the bottom surface should be laminar
lamLengthBottom = 0.7 #position where MPR starts at bottom -> lower end of bucket
#Ramp Bottom Side against lam seperation
blendingLamBottom = 0.05 #reach of ramp into lam area -> increase when laminar separation occurs at MPR
blendingMPRBottom = 0.01 #ramp into MPR -> to combat laminar separation 
#MPR bottom side 
shapeMPRBottom = 0.23
strengthMPRBottom = 0.7

##Top Side
lamTop = 0.9 #Target Ca for laminar flow at Top -> upper end of bucket
lamLengthTop = 0.6 #laminar length
#Ramp Top Side against lam seperation
blendingLamTop = 0.05 #reach of ramp into lam area -> increase when laminar separation occurs at MPR
blendingMPRTop = 0.01 #ramp into MPR -> to combat laminar separation 
#MPR Top side 
shapeMPRTop = -1
strengthMPRTop = 0.7
#alpah* increase against suction peak -> increases alpha* smothly (according to a polynom) towards the LE to combat suction peaks
alphaLE = 14 #the alpha* we will have at the LE
startIncrease = 0.2  #the c at which the increase starts
numberPoints = 8    #amount of points to use
alphaMax = 18 #maximum alpha we will reach between Increase Start and LE
alphaMaxLocation = 0.08

adaptSide = "upper" #which side to iterate for MPR choices are "upper", "lower", "mixed", "none"


##GenerateFile
foil1 = EpplerFoil.AirFoil("entwurf.dat","TF",421)
foil1.BeginnFile()

#generate upper side alpha*
#fit function to combat LE
#func is a*x**3 + b*x**2 + c*x + d + e*x**4 + g*x**5

opti = asb.Opti() #erstellt optimizer environment

a = 0
b = 0
c = 0
d = opti.variable(init_guess=1)
e = opti.variable(init_guess=1)
g = opti.variable(init_guess=1)

#function to interpolate:
def f(x,a,b,c,d,e,g):
    return a + b*x + c*x**2 + d*x**3 + e*x**4 + g*x**5
def df(x,a,b,c,d,e,g):
    return b + 2*c**x + 3*d*x**2 + 4*e*x**3 + 5*g*x**4
def ddf(x,a,b,c,d,e,g):
    return 2*c + 6*d*x + 12*e*x**2 + 20*g*x**3

#transformation
delta = startIncrease
xMax = (startIncrease - alphaMaxLocation) / delta
xLE = 1
aMax = alphaMax - (lamTop/0.11)
aLE = alphaLE - (lamTop/0.11)

opti.subject_to(f(xMax,a,b,c,d,e,g) == aMax)
opti.subject_to(df(xMax,a,b,c,d,e,g) == 0)
opti.subject_to(f(xLE,a,b,c,d,e,g) == aLE)

sol = opti.solve()
d = sol.value(d)
e = sol.value(e)
g = sol.value(g)


y = np.zeros(numberPoints + 1)
x = np.zeros(numberPoints + 1)
for i in range(0,numberPoints):
    x[i] = startIncrease - delta * i / numberPoints
    y[i] = f(i/numberPoints,a,b,c,d,e,g) + (lamTop/0.11)
x[-1] = 0
y[-1] = alphaLE

foil1.WriteUpperC_Alpha(x,y)
foil1.MPRUpper(lamLengthTop,0.98,shapeMPRTop,strengthMPRTop)
foil1.RampUpper(blendingLamTop,blendingMPRTop)

#Bottom Side
foil1.WriteLowerC_Alpha([1],[(lamBottom/0.11)])
foil1.MPRLower(lamLengthBottom,0.98,shapeMPRBottom,strengthMPRBottom)
foil1.RampLower(blendingLamBottom,blendingMPRBottom)
foil1.writeMPR(adaptSide,0.7,0)


#Analysis
foil1.inviscidCalc_byIncrements(0,3,7)
foil1.naturalTransitionCalc(3000)
foil1.visousCalc(0,15,15)
foil1.CloseFile()
foil1.ExecuteEppler()





