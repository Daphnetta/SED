import matplotlib.pyplot as plt
from math import log, sqrt, sin, cos, asin, radians, degrees, pi, exp
from scipy import integrate
import numpy as np

M_star = 2*10**33
R_star = 2*7*10**10
d = 150* 3.086 * 10**8
delta = (10**(-2))*R_star
M_dot = M_star*10**(-8) 
i = 1
T_eff_star = 4000
G = 6.67*10**(-8)
sigma = 5.67*10**(-5)
R_D = 3*R_star

m_dot = 1
m = 1
h = 6.626* 10**(-27)
c = 3*10**10
k_B = 1.38*10**(-16)

lst = []
ast = []
bst = []

def Disk():
	for lam in np.arange (3*10**(-7) , 10**(-2), 10**(-6)) :
		print(lam)
		lst.append(log(lam))
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

	return bst, lst

def Star():
	for lam in np.arange(10**(-6) , 10**(-2), 3*10**(-7)):
		print(lam)
		lst.append(log(lam))
		f_all_1 = (0,0)
		for r in np.arange(0, R_star, R_star/100):

			T_eff = lambda r: 150*(r/(1.5*10**13))**-0.75
			#print(T_eff)
			f = lambda r: ((2*pi*h*c**2)/(lam**5))*(1/(exp((h*c)/(lam*k_B*T_eff(r)))-1))
			r_1 = integrate.quad(f,  R_star, R_D)
			f_all_1 = (f_all_1[0]+r_1[0], sqrt(f_all_1[1]**2+r_1[1]**2))
			#print(f_all_1[0])
			#print(f_all_1)
			#ast.append(f_all[0])		
			#print(f_all_1[0])
		bst.append(log(lam*(f_all_1[0]))) 
	return bst, lst


ax = plt.gca()
Star()
plt.plot(lst, bst,':m', label='Star')
lst = []
ast = []
bst = []
Disk()
plt.plot(lst, bst,'--g', label='Disk')


ax.set_ylim([0, 65])
ax.set_xlim([-15, -5])
plt.xlabel('$\\log \\lambda\; [ \mathrm{cm}$]')
plt.ylabel('$\\log \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.title('Results')
plt.legend()

#Star()
#ax = plt.gca()
#plt.plot(lst,bst)

#ax.set_ylim([35, 65])
plt.show()
