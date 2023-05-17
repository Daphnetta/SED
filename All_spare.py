import matplotlib.pyplot as plt
from math import log, log10, sqrt, sin, cos, asin, radians, degrees, pi, exp
from scipy import integrate
import numpy as np
import sys
sys.path.append("./")
from disk_r_structure import Disk
from const import pc, M_sun, au, k, R_sun, year, h, c
m = 1
mdot = 1
alpha = 1
l_star = 2
r_au = 1
disk = Disk(m, mdot, alpha, l_star)
disk.calc_const()

M_star = M_sun            # g
R_star = 2*R_sun          # cm
d = 150.0*pc              # cm
M_dot = 1.0e-8*M_sun/year # g/c 
R_D = 100.0*au            # cm
i = radians(0)            # grad
T_eff_star = 4e3          # K

lambda_min = 1e-6         # cm
lambda_max = 1e-1         # cm
N = 200
log_lambdas = np.linspace(log10(lambda_min), log10(lambda_max), N)
lambdas = 10**log_lambdas

x_min = 3
x_max = R_D / R_star

def Teff(z, r_au):
	if z == 1:
		return disk.Teff_v(r_au)
	elif z == 2:
		return disk.Teff_irr(r_au)

def gamma_0(x, i):
    if ((1 < x) and (x < 1.0 / cos(i))):
        return asin(sqrt(1 - x**(-2)) / sin(i))
    else:
        return 0.5 * pi

def subint(x, lam):
	arg = h*c / (lam*k*Teff(z, r_au))
	if (arg < 700):
		return x*(pi + 2.0 * gamma_0(x, i)) / (exp(arg) - 1)
	else:
		return x*(pi + 2.0 * gamma_0(x, i)) / (arg)

def Disk1(z):
	lamst = []	# пустой список для значений lambda
	lamfst = [] # пустой список для значений lambda*F_lambda
	for lam in lambdas:
		lamst.append(lam)
		#print(lam)
		b = ((2*h*c**2)/(lam**5))*((R_star/d)**2)* cos(i) # константа к I_lambda
		f_all_1 = (0,0)
		Interg = lambda x: b*subint(x, lam)
		r_1 = integrate.quad(Interg, x_min, x_max)
		f_all_1 = (f_all_1[0]+r_1[0], sqrt(f_all_1[1]**2+r_1[1]**2))	
		#print(rs)
		if lam < 3e-4:
			lamfst.append(0)
		else:
			lamfst.append(lam*(f_all_1[0])) 

	#plt.plot(lamst, lamfst,'--g', label='Disk')
	return lamfst

def Star():
	lamst = [] 	# пустой список для значений lambda
	lamfst = [] # пустой список для значений lambda*F_lambda
	for lam in lambdas:
		lamst.append(lam) # запись в список lambda
		b = (2*pi*h*c**2)/(lam**5) # константа к I_lambda
		с = (pi/2)*(1 + cos(i))*(R_star/d)**2 # константа к F_lambda
		I = b*(1/(exp((h*c)/(lam*k*T_eff_star))-1))
		f = с*I
		lamfst.append(lam*f) # запись в список lambda*F_lambda
	return lamfst

ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim([10**-15, 10**-4])
#ax.set_xlim([10**-6, 10**-1])
lamst = []
for lam in lambdas:
	lamst.append(lam)
lamfst = Star()
plt.plot(lamst, lamfst,'b', label='Star')
z=1
lamfall_1 = np.array([Star(),Disk1(z)]).sum(axis=0).tolist()
z=2
lamfall_2 = np.array([Star(),Disk1(z)]).sum(axis=0).tolist()
#print(lamfall)
plt.plot(lamst, lamfall_1,':g', label='T_eff_v')
plt.plot(lamst, lamfall_2,':m', label='T_eff_irr')
plt.xlabel('$\\log \\lambda\; [ \mathrm{cm}$]')
plt.ylabel('$\\log \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.title('r_au=10')
plt.legend()
plt.show()

