import matplotlib.pyplot as plt
from math import sqrt, sin, cos, asin, radians, degrees, pi, exp, log
from scipy import integrate
import numpy as np
#print(sin(radians(90)))

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

lst = []
ast = []

bst = []
for lam in np.arange (10**(-6) , 10**(-3), 10**(-6)) :
	lst.append(log(lam))
	#accretion disk
	f_all_1 = (0,0)
	f_all_2 = (0,0)
	for r in np.arange (R_star, R_D, R_star/100) :
		Phi = lambda r: asin(R_star/r)
		F_V = lambda r: ((3*G*M_star*M_dot)/(8*pi*r**3))*(1-sqrt(R_star/r))
		F_A = lambda r: ((sigma/pi)*(T_eff_star)**4)*(Phi(r)-(sin(2*Phi(r)))/2)
		T_D = lambda r: ((F_V(r) + F_A(r))/sigma)**0.25
		x = r/R_star
		if r>R_star and r<((1/cos(radians(i)))*R_star):
			gamma_0 = lambda x: asin(((1-(x)**(-2))**0.5)/sin(radians(i)))
			Interg = lambda x: ((2*h*c**2)/(lam**5))*((R_star/d)**2)* cos(radians(i)) * ((pi+2*gamma_0(x))/(exp((h*c)/(lam*k_B*T_D(r/R_star)))-1))*x
			r_1 = integrate.quad(Interg, 1, R_D/R_star)
			f_all_1 = (f_all_1[0]+r_1[0], sqrt(f_all_1[1]**2+r_1[1]**2))	
		elif r>((1/cos(radians(i)))*R_star) :
			gamma_0 = pi/2
			Interg = lambda x: ((2*h*c**2)/(lam**5))*((R_star/d)**2)* cos(radians(i)) * ((pi+2*gamma_0)/(exp((h*c)/(lam*k_B*T_D(r)))-1))*x
			r_2 = integrate.quad(Interg, 1, R_D/R_star)
			f_all_2 = (f_all_2[0]+r_2[0], sqrt(f_all_2[1]**2+r_2[1]**2))
	bst.append(log(lam*(f_all_1[0]+f_all_2[0]))) 

print(lst)
print(bst)
ax = plt.gca()
plt.plot(lst,bst)

ax.set_ylim([35, 65])
plt.show()
