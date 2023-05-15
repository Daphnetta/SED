import matplotlib.pyplot as plt
from math import log, log10, sqrt, sin, cos, asin, radians, degrees, pi, exp
from scipy import integrate
import numpy as np

M_star = 2e33      # g
R_star = 1.4e11    # cm
d = 4.629e10       # cm
M_dot = M_star*1e-8  # g/c 
R_D = 2*R_star     # cm
m_dot = 1
m = 1
i = 1              # grad
T_eff_star = 4e3   # K

G = 6.67e-8        # cm3*g-1*s-2
sigma = 5.67e-5    # g*cm-3**K-4
h = 6.626e-27      # g*cm2*s-1
c = 3e10           # cm/s
k_B = 1.38e-16     # g*cm2*s-2*K-1

lambda_min = 1e-6  # cm
lambda_max = 1e-3  # cm
N = 200
log_lambdas = np.linspace(log10(lambda_min), log10(lambda_max), N)
lambdas = 10**log_lambdas

r_min = 0  	        # cm
r_max = R_star      # cm
N_r = 1e8
r_s = np.arange(r_min, r_max, N_r)

for lam in lambdas:
	print(lam) # для отслеживания процесса цикла
	lamst.append(lam)
	b = (2*pi*h*c**2)/(lam**5) # константа к I_lambda
	с = (pi/2)*(1 + cos(radians(i)))*(R_star/d)**2 # константа к F_lambda
	I = b*(1/(exp((h*c)/(lam*k_B*T_eff_star))-1))
	f = с*I
	lamfst.append(lam*f) 

print(lamst)
print(lamfst)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
plt.plot(lamst, lamfst,':m', label='Star') # построение графика
plt.xlabel('$\\log \\lambda\; [ \mathrm{cm}$]')
plt.ylabel('$\\log \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.title('Central star')
plt.legend()
plt.show()
