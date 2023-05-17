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
R_D = 1.5e15            # cm
m_dot = 1
m = 1
i = radians(0)                     # grad
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

lamst = [] 	              # пустой список для значений lambda
lamfst = []               # пустой список для значений lambda*F_lambda

for lam in lambdas:
	lamst.append(lam) # запись в список lambda
	b = (2*pi*h*c**2)/(lam**5) # константа к I_lambda
	с = (pi/2)*(1 + cos(i))*(R_star/d)**2 # константа к F_lambda
	I = b*(1/(exp((h*c)/(lam*k_B*T_eff_star))-1))
	f = с*I
	lamfst.append(lam*f) # запись в список lambda*F_lambda

ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim([10**-15, 10**-4])
plt.plot(lamst, lamfst,':m', label='Star') # построение графика
plt.xlabel('$\\log \\lambda\; [ \mathrm{cm}$]')
plt.ylabel('$\\log \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.title('Central star')
plt.legend()
plt.show()
