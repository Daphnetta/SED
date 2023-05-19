import matplotlib.pyplot as plt
from math import log, log10, sqrt, sin, cos, asin, radians, pi, exp
from scipy import integrate
import numpy as np
import sys
sys.path.append("./")
from disk_r_structure import Disk
from const import pc, M_sun, au, k, R_sun, h, c
m = 1
mdot = 1  # = 1 for _md8_; or = 0.1 for _md9_ in datafile name
alpha = 0.01 # = 1 for alpha001; = 0.1 for alpha0001; = 0.01 for alpha00001 in datafile name
l_star = 2
disk = Disk(m, mdot, alpha, l_star)
disk.calc_const()
disk.import_data("./disk_data/alpha00001_md8_amft0001_cr17_xr30/solution.dat")

M_star = M_sun            # g
R_star = 2*R_sun          # cm
d = 150.0*pc              # cm
R_D = 100.0*au            # cm
i = radians(0)            # grad
T_eff_star = 4e3          # K

lambda_max = 1e-1         # cm
N = 200

def make_lambda_grid(l_min,l_max, N):
	log_lambdas = np.linspace(log10(l_min), log10(l_max), N)
	lambdas = 10**log_lambdas
	return lambdas

print("r_in  = ", disk.Get_r_in(), "au")
print("r_out = ", R_D / au, "au")

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

def Disk1(lambda_min, xmin, xmax, z):
	lamst = []	# пустой список для значений lambda
	lamfst = [] # пустой список для значений lambda*F_lambda
	lambdas = make_lambda_grid(lambda_min, lambda_max, N)

	for lam in lambdas:
		lamst.append(lam)
		#print(lam)
		b = ((2*h*c**2)/(lam**5))*((R_star/d)**2)* cos(i) # константа к I_lambda
		f_all_1 = (0,0)
		Interg = lambda x: b*subint(z, x, lam)
		r_1 = integrate.quad(Interg, xmin, xmax)
		f_all_1 = (f_all_1[0]+r_1[0], sqrt(f_all_1[1]**2+r_1[1]**2))	
		lamfst.append(lam*(f_all_1[0])) 

	#plt.plot(lamst, lamfst,'--g', label='Disk')
	#s_log = list(zip([log(w) for w in lamst], [log(w) for w in lamfst]))
	#k_xy_log = list([(t[1][1] - t[0][1]) / (t[1][0] - t[0][0]) for t in zip(s_log, s_log[1:])])
	#print(k_xy_log[155:185])
	#plt.plot(lamst[155:185], lamfst[155:185], 'b')
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

	#s_log = list(zip([log(w) for w in lamst], [log(w) for w in lamfst]))
	#k_xy_log = list([(t[1][1] - t[0][1]) / (t[1][0] - t[0][0]) for t in zip(s_log, s_log[1:])])
	#print(k_xy_log[130:151])
	#plt.plot(lamst[130:151], lamfst[130:151], 'b')
	return [lamst, lamfst]

#ax = plt.gca()
#ax.set_xscale("log")
#ax.set_yscale("log")
#ax.set_ylim([10**-15, 10**-4])
#ax.set_xlim([10**-6, 10**-1])

[lamst, lamfst] = Star(1e-6)
#plt.plot(lamst, lamfst,':y', label='Star')

x_min = 5*au/R_star
x_max = R_D / R_star
disk_data_1 = Disk1(5e-4,  x_min, x_max, 1)
disk_data_2 = Disk1(1e-4,  x_min, x_max, 2)

x_min = disk.Get_r_in() * au / R_star
x_max = R_D / R_star
disk_data_3 = Disk1(6e-4, x_min, x_max, 3)

#plt.plot(disk_data_1[0], disk_data_1[1],':g', label='T_eff_v')
#plt.plot(disk_data_2[0], disk_data_2[1],':m', label='T_eff_irr')
#plt.plot(disk_data_3[0], disk_data_3[1],'--r', label='T_eff_num')

#plt.xlabel('$ \\lambda\; [ \mathrm{cm}$]')
#plt.ylabel('$ \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
#plt.legend()
#plt.show()


#All = np.array([lamst, lamfst])
#All1 = np.append(All, disk_data_1,axis = 0)
#All2 = np.append(All1, disk_data_2,axis = 0)
#All3 = np.append(All2, disk_data_3,axis = 0)
#np.savetxt("./disk_data/alpha00001_md9_amft0001_cr17_xr30/mgd_result.dat",np.transpose(All3), delimiter=' ', newline='\n', header='1:lamst 2:lamfst 3:disk_data_1[0] 4:disk_data_1[1] 5:disk_data_1[0] 6:disk_data_1[1] 7:disk_data_3[0] 8:disk_data_3[1] ') 

def Result():

	data = np.loadtxt("./disk_data/alpha00001_md8_amft0001_cr17_xr30/result.dat", skiprows=1)
	lamst = data[:, 0]
	lamfst = data[:, 1]
	disk_data_1 = [[],[]]
	disk_data_2 = [[],[]]
	disk_data_3 = [[],[]]
	disk_data_4 = [[],[]]
	disk_data_1[0] = data[:, 2]
	disk_data_1[1] = data[:, 3]
	disk_data_2[0] = data[:, 4]
	disk_data_2[1] = data[:, 5]
	disk_data_3[0] = data[:, 6]
	disk_data_3[1] = data[:, 7]
	ax = plt.gca()
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_ylim([10**-15, 10**-4])
	ax.set_xlim([10**-6, 10**-1])

	plt.plot(lamst, lamfst,'y', label='Star')
	plt.plot(disk_data_1[0], disk_data_1[1],'g', label='T_eff_v')
	plt.plot(disk_data_2[0], disk_data_2[1],'m', label='T_eff_irr')
	plt.plot(disk_data_3[0], disk_data_3[1],'r', label='T_eff_num')

	data = np.loadtxt("./disk_data/alpha00001_md8_amft0001_cr17_xr30/mgd_result.dat", skiprows=1)
	disk_data_4[1] = data[:, 7]
	plt.plot(disk_data_3[0], disk_data_4[1],':r', label='T_eff_num_mgd')

	plt.xlabel('$ \\lambda\; [ \mathrm{cm}$]')
	plt.ylabel('$ \\lambda F_\\lambda \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
	plt.legend()
	plt.show()


Result()
