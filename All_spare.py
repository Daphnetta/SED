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
disk = Disk(m, mdot, alpha, l_star)
disk.calc_const()
disk.import_data("solution.dat")

M_star = M_sun            # g
R_star = 2*R_sun          # cm
d = 150.0*pc              # cm
M_dot = 1.0e-8*M_sun/year # g/c 
R_D = 100.0*au            # cm
i = radians(0)            # grad
T_eff_star = 4e3          # K

lambda_max = 1e-1         # cm
N = 200

def make_lambda_grid(l_min,l_max, N):
	log_lambdas = np.linspace(log10(l_min), log10(l_max), N)
	lambdas = 10**log_lambdas
	return lambdas

x_min = 1 * au / R_star
x_max = R_D / R_star

def Teff(z, r_au):
	if z == 1:
		return disk.Teff_v(r_au)
	elif z == 2:
		return disk.Teff_irr(r_au)
	elif z == 3:
		return disk.Tnum(r_au)
	else:
		return disk.Teff_irr(r_au)

def gamma_0(x, i):
	if ((1 < x) and (x < 1.0 / cos(i))):
		return asin(sqrt(1 - x**(-2)) / sin(i))
	else:
		return 0.5 * pi

def subint(z, x, lam):
	arg = h * c / (lam * k * Teff(z, x * R_star / au))
	if (arg < 700):
		return x * (pi + 2.0 * gamma_0(x, i)) / (exp(arg) - 1)
	else:
		return x * (pi + 2.0 * gamma_0(x, i)) / (arg)

def Disk1(lambda_min, z):
	lamst = []	# пустой список для значений lambda
	lamfst = [] # пустой список для значений lambda*F_lambda
	lambdas = make_lambda_grid(lambda_min, lambda_max, N)

	for lam in lambdas:
		lamst.append(lam)
		#print(lam)
		b = ((2*h*c**2)/(lam**5))*((R_star/d)**2)* cos(i) # константа к I_lambda
		f_all_1 = (0,0)
		Interg = lambda x: b*subint(z, x, lam)
		r_1 = integrate.quad(Interg, x_min, x_max)
		f_all_1 = (f_all_1[0]+r_1[0], sqrt(f_all_1[1]**2+r_1[1]**2))	
		#print(rs)
		if lam < 3e-5:
			lamfst.append(0)
		else:
			lamfst.append(lam*(f_all_1[0])) 

	#plt.plot(lamst, lamfst,'--g', label='Disk')
	return [lamst, lamfst]

def Star(lambda_min):
	lamst = [] 	# пустой список для значений lambda
	lamfst = [] # пустой список для значений lambda*F_lambda
	lambdas = make_lambda_grid(lambda_min, lambda_max, N)
	for lam in lambdas:
		lamst.append(lam) # запись в список lambda
		b = (2*pi*h*c**2)/(lam**5) # константа к I_lambda
		с = (pi/2)*(1 + cos(i))*(R_star/d)**2 # константа к F_lambda
		I = b*(1/(exp((h*c)/(lam*k*T_eff_star))-1))
		f = с*I
		lamfst.append(lam*f) # запись в список lambda*F_lambda
	return [lamst, lamfst]

ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim([10**-15, 10**-4])
#ax.set_xlim([10**-6, 10**-1])

[lamst, lamfst] = Star(1e-6)
plt.plot(lamst, lamfst,'b', label='Star')

disk_data_1 = Disk1(5e-4, 1)
disk_data_2 = Disk1(1e-4, 2)
disk_data_3 = Disk1(1e-3, 3)

plt.plot(disk_data_1[0], disk_data_1[1],':g', label='T_eff_v')
plt.plot(disk_data_2[0], disk_data_2[1],':m', label='T_eff_irr')
plt.plot(disk_data_3[0], disk_data_3[1],'--r', label='T_eff_num')

plt.xlabel('$\\log \\lambda\; [ \mathrm{cm}$]')
plt.ylabel('$\\log \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.title('r_au=10')
plt.legend()
plt.show()

