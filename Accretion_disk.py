import matplotlib.pyplot as plt
from math import log, log10, sqrt, sin, cos, asin, radians, degrees, pi, exp
from scipy import integrate
import numpy as np

M_sun = 1.99e33           # g
R_sun = 6.96e10           # cm 
M_star = M_sun            # g
R_star = 2*R_sun          # cm
d = 4.629e20              # cm
M_dot = M_star*3.1536e-1  # g/c 
R_D = 1.5e15              # cm
m_dot = 1
m = 1
i = radians(0)                     # grad
T_eff_star = 4e3          # K

G = 6.67e-8               # cm3*g-1*s-2
sigma = 5.67e-5           # g*cm-3**K-4
h = 6.626e-27             # g*cm2*s-1
c = 3e10                  # cm/s
k_B = 1.38e-16            # g*cm2*s-2*K-1

lambda_min = 1e-12         # cm
lambda_max = 1e-1         # cm
N = 200
log_lambdas = np.linspace(log10(lambda_min), log10(lambda_max), N)
lambdas = 10**log_lambdas

x_min = 1
x_max = R_D / R_star

lamst = [] 	              # пустой список для значений lambda
lamfst = []               # пустой список для значений lambda*F_lambda

def T_D_eff(x):
	return 150*(x/(1.5*10**13))**-0.75

def gamma_0(x, i):
    if ((1 < x) and (x < 1.0 / cos(i))):
        return asin(sqrt(1 - x**(-2)) / sin(i))
    else:
        return 0.5 * pi

def subint(x, lam):
    arg = h*c / (lam*k_B*T_D_eff(x*R_star))
    if (arg < 700):
        return x*(pi + 2.0 * gamma_0(x, i)) / (exp(arg) - 1)
    else:
        return x*(pi + 2.0 * gamma_0(x, i)) / (arg)

for lam in lambdas:
	lamst.append(lam)
	print(lam)
	b = ((2*h*c**2)/(lam**5))*((R_star/d)**2)* cos(i) # константа к I_lambda
	f_all_1 = (0,0)
	Interg = lambda x: b*subint(x, lam)
	r_1 = integrate.quad(Interg, x_min, x_max)
	f_all_1 = (f_all_1[0]+r_1[0], sqrt(f_all_1[1]**2+r_1[1]**2))	
	#print(rs)
	lamfst.append(lam*(f_all_1[0])) 

ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_ylim([10**7, 10**16])
plt.plot(lamst, lamfst,'--g', label='Disk') # построение графика
plt.xlabel('$\\log \\lambda\; [ \mathrm{cm}$]')
plt.ylabel('$\\log \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.title('Accretion disk')
plt.legend()
plt.show()