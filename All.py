import matplotlib.pyplot as plt
from math import log, log10, sqrt, sin, cos, asin, radians, degrees, pi, exp
from scipy import integrate
import numpy as np

M_sun = 1.99e33           # g
R_sun = 6.96e10           # cm 
M_star = M_sun            # g
R_star = 2*R_sun          # cm
d = 4.629e10              # cm
M_dot = M_star*3.1536e-1  # g/c 
R_D = 2*R_star            # cm
m_dot = 1
m = 1
i = 1                     # grad
T_eff_star = 4e3          # K

G = 6.67e-8               # cm3*g-1*s-2
sigma = 5.67e-5           # g*cm-3**K-4
h = 6.626e-27             # g*cm2*s-1
c = 3e10                  # cm/s
k_B = 1.38e-16            # g*cm2*s-2*K-1

lambda_min = 1e-6         # cm
lambda_max = 1e-1         # cm
N = 200
log_lambdas = np.linspace(log10(lambda_min), log10(lambda_max), N)
lambdas = 10**log_lambdas

r_min = R_star  	      # cm
r_max = R_D               # cm
N_r = 1e8
rs = (r_min, r_max, N_r)

def T_D_eff(r):
	return 150*(r/(1.5*10**13))**-0.75

def Disk():
	lamst = [] 	# пустой список для значений lambda
	lamfst = [] # пустой список для значений lambda*F_lambda
	for lam in lambdas:
		lamst.append(lam)
		b = ((2*h*c**2)/(lam**5))*((R_star/d)**2)* cos(radians(i)) # константа к I_lambda
		f_all_1 = (0,0)
		f_all_2 = (0,0)
		for r in rs :
			Phi = lambda r: asin(R_star/r)
			x = r/R_star
			if r>R_star and r<((1/cos(radians(i)))*R_star):
				gamma_0 = lambda x: asin(((1-(x)**(-2))**0.5)/sin(radians(i)))
				Interg = lambda x: b*((pi+2*gamma_0(x))/(exp((h*c)/(lam*k_B*T_D_eff(r/R_star)))-1))*x
				r_1 = integrate.quad(Interg, 1, R_D/R_star)
				f_all_1 = (f_all_1[0]+r_1[0], sqrt(f_all_1[1]**2+r_1[1]**2))	
			elif r>((1/cos(radians(i)))*R_star) :
				gamma_0 = pi/2
				Interg = lambda x: b*((pi+2*gamma_0)/(exp((h*c)/(lam*k_B*T_D_eff(r)))-1))*x
				r_2 = integrate.quad(Interg, 1, R_D/R_star)
				f_all_2 = (f_all_2[0]+r_2[0], sqrt(f_all_2[1]**2+r_2[1]**2))
		lamfst.append(lam*(f_all_1[0]+f_all_2[0])) 
	plt.plot(lamst, lamfst,'--g', label='Disk')
	return lamst, lamfst

def Star():
	lamst = [] 	# пустой список для значений lambda
	lamfst = [] # пустой список для значений lambda*F_lambda
	for lam in lambdas:
		lamst.append(lam) # запись в список lambda
		b = (2*pi*h*c**2)/(lam**5) # константа к I_lambda
		с = (pi/2)*(1 + cos(radians(i)))*(R_star/d)**2 # константа к F_lambda
		I = b*(1/(exp((h*c)/(lam*k_B*T_eff_star))-1))
		f = с*I
		lamfst.append(lam*f) # запись в список lambda*F_lambda
	plt.plot(lamst, lamfst,':m', label='Star')
	return lamst, lamfst

ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim([10**7, 10**12])
Star()
Disk()
plt.xlabel('$\\log \\lambda\; [ \mathrm{cm}$]')
plt.ylabel('$\\log \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.title('')
plt.legend()
plt.show()

