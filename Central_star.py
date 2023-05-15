import matplotlib.pyplot as plt
from math import log, log10, sqrt, sin, cos, asin, radians, degrees, pi, exp
from scipy import integrate
import numpy as np

M_star = 2e33      # g
R_star = 1.4e11    # cm
d = 4.629e10       # cm
M_dot = M_star/1e8  # g/c 
R_D = 2*R_star     # cm
m_dot = 1
m = 1
i = 1              # grad
T_eff_star = 4e3   # K

G = 6.67/1e8        # cm3*g-1*s-2
sigma = 5.67/1e5    # g*cm-3**K-4
h = 6.626/1e27      # g*cm2*s-1
c = 3e10           # cm/s
k_B = 1.38/1e16     # g*cm2*s-2*K-1
#f_all = 0.0
lamst = [] 	# пустой список для значений lambda
lamfst = [] # пустой список для значений lambda*F_lambda
#f_all_1 = []
#f_all_2 = []
f_all_1 = [0,0] # список для значений вычисленных интегралов

lambda_min = 1/1e6  # cm
lambda_max = 1/1e3  # cm
N = 200
log_lambdas = np.linspace(log10(lambda_min), log10(lambda_max), N)
lambdas = 10**log_lambdas

r_min = 0  	        # cm
r_max = R_star      # cm
N_r = 1e8
r_s = np.arange(r_min, r_max, N_r)

for lam in lambdas:
	print(lam) # для отслеживания процесса цикла
	lamst.append(log10(lam)) # запись в список log10(lambda)
	f_all_1 = (0,0) # кортеж
	b = (2*pi*h*c**2)/(lam**5) # константа к I_lambda(r)
	с = (pi/2)*(1 + cos(radians(i)))*(R_star/d)**2 # константа к F_lambda(r)
	for r in r_s:
		T_eff = lambda r: 150*(r/(1.5e13))**-0.75 # описание функции температуры T_eff(r)
		#print(T_eff)
		I = lambda r: b*(1/(exp((h*c)/(lam*k_B*T_eff(r)))-1)) # описание функции I_lambda(r)
		f = lambda r: с*I(r) # описание функции F_lambda(r)
		r_1 = integrate.quad(f,  0, R_star) # интегрирование функции F_lambda(r) в пределах от 0 до R_star
		f_all_1 = ((f_all_1[0]+r_1[0]), (sqrt(f_all_1[1]**2+r_1[1]**2))) # запись результатов интегрирования: 1 - значение интеграла, 2 - ошибка
	lamfst.append(log10(lam*(f_all_1[0]))) # запись в список log(lambda*F_lambda) 

print(lamst)
print(lamfst)
ax = plt.gca()
plt.plot(lamst, lamfst,':m', label='Star') # построение графика
plt.xlabel('$\\log \\lambda\; [ \mathrm{cm}$]')
plt.ylabel('$\\log \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
plt.title('Central star')
plt.legend()
plt.show()
