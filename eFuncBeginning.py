import aerosandbox as asb
import aerosandbox.numpy as np

n = 5 #amount of points to give
i0 = 26.1 #starting point
i1 = 31.1 #endpoint
alpha0 = 8.3 #alpha at start of increase
alpha1 = 20 #alpha at the end of increase

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

opti.subject_to(f(i0,a,b,c,d) == alpha0)
opti.subject_to(f(i1,a,b,c,d) == alpha1)
opti.subject_to(df(i0,a,b,c,d) == 0)
opti.subject_to(ddf(i0,a,b,c,d) == 0)

sol = opti.solve()
a = sol.value(a)
b = sol.value(b)
c = sol.value(c)
d = sol.value(d)

xarr = np.linspace(i0,i1,num=n)
for i in range(0,n):
    x = xarr[i]
    y = a*x**3 + b*x**2 + c*x + d
    print(str(round(x,1)) + " " + str(round(y,1)), end= ' ')
