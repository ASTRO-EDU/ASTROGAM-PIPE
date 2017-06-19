from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import scipy.optimize
import scipy.stats


#THETA E PHI TEORICI
theta_th = 30.
phi_th = 225.

energy = 100.  #MeV
#./theta'+str(int(theta_th))+'phi'+str(int(phi_th))+'onlycal/kalman_output/eAST211112021_CLUSTER_100000ph_Point_MONO_UNI_100MeV_'+str(int(theta_th))+'_'+str(int(phi_th))+'_PFIND_1_5.dat
#eAST211112021_CLUSTER_100000ph_Point_MONO_UNI_100MeV_30_225.0_PFIND_1_5.dat
with open('./theta'+str(int(theta_th))+'phi'+str(int(phi_th))+'onlycal/kalman_output/eAST211112021_CLUSTER_100000ph_Point_MONO_UNI_100MeV_'+str(int(theta_th))+'_'+str(int(phi_th))+'_PFIND_1_5.dat') as f:
	f=[x.strip() for x in f if x.strip()]
	data=[tuple(map(float,x.split())) for x in f[0:]]
	event=[x[0] for x in data]
	theta=[x[1] for x in data]
	phi=[x[2] for x in data]
	
sph_dist_list = []	
i = 0
while i< len(event)-1:
	#print('ev',int(event[i]), 'th',theta[i], 'ph',phi[i])
	if ((theta[i] == theta_th) and (phi[i] == phi_th)):
		print('distanza sferica = 0')
	if ((str(theta[i]) == 'nan') or (str(phi[i]) == 'nan')):
		i = i +1
	if ((str(theta[i]) == 'nan') and (str(phi[i]) == 'nan')):
		i = i +1

	#DISTANZA SFERICA
	sph_dist = math.acos((math.sin(theta[i])*math.sin(theta_th)*math.cos(phi[i] - phi_th) + math.cos(theta[i])*math.cos(theta_th)))
	#print(i, sph_dist)
	sph_dist_list.append(sph_dist)

	i = i +1

#ORDINAMENTO VETTORE
sph_dist_list = np.sort(sph_dist_list)

#68% ----> 1 sigma
j = int(len(sph_dist_list)/100.*68.)
#print('jesimo elemento al 68 percento', j)
psf = sph_dist_list[j]
errore = int(math.sqrt(j))
#print('errore sulla j',errore)
err_inf = j-errore
err_sup = j+errore
#print('errori inf e sup sulle j', err_inf, err_sup)
psf_inf = sph_dist_list[err_inf]
psf_sup = sph_dist_list[err_sup]
#print('corrispondenti psf', psf_inf, psf_sup)
psf_error_inf = psf-psf_inf
psf_error_sup = -(psf-psf_sup)
#print(psf_error_inf, psf_error_sup)
media_errore = (psf_error_inf+psf_error_sup)/2.
#print('errore medio',media_errore)

#CONTROLLO VETTORE ORDINATO
#print(len(sph_dist_list))
#i = 0
#while i< len(sph_dist_list):
#	print(i, sph_dist_list[i])
#	i = i +1
	
print('La PSF dello strumento e'' '+str(psf)+' +- '+str(media_errore))+' gradi '	

#j = 0
#while j < len(sph_dist_list):
#	counts = conteggi/(math.pi*(sph_dist_list[i])**2)

n_bins = 20
# creare l istogramma e normalizzare
N_array, bin_array = np.histogram(sph_dist_list, bins=n_bins) 
norm_N = np.zeros(len(N_array))                                          # creo il vettore vuoto per gli N normalizzati
err_norm_N = np.zeros(len(N_array))                                # creo il vettore vuoto per gli errori propagati in asse Y
err_theta = np.zeros(len(N_array))                               # creo il vettore vuoto per gli errori in asse X
theta_array = np.zeros(len(N_array))                          # creo il vettore vuoto per l array di theta in asse X
y_king_arr = np.zeros(len(sph_dist_list))


for jn in xrange(len(N_array)):
	err_theta[jn] = (bin_array[jn+1] - bin_array[jn])/2.
	theta_array[jn] = bin_array[jn] + err_theta[jn]
	theta_bin = bin_array[jn+1] - bin_array[jn]  # la larghezza del bin
	area = math.pi*((sph_dist_list[jn])**2)
	theta_bin_norm = theta_bin*area
	norm_N[jn] = float(N_array[jn])/theta_bin_norm  # divido per la larghezza del bin
	err_norm_N[jn] = (np.sqrt(float(N_array[jn])))/theta_bin_norm  # stessa cosa per l errore
	


norm_N_bis = []
err_norm_N_bis = []
err_theta_bis = []
theta_array_bis = []
theta_array_bis_x_arr = []
i = 0
while i < len(N_array)-1:
	norm_N_bis.append(norm_N[i])
	err_norm_N_bis.append(err_norm_N[i])
	err_theta_bis.append(err_theta[i])
	theta_array_bis.append(theta_array[i])
	i = i +1
	
norm_N_bis = np.array(norm_N_bis)
err_norm_N_bis = np.array(err_norm_N_bis)
err_theta_bis = np.array(err_theta_bis)
theta_array_bis = np.array(theta_array_bis)

s_p_e = math.sqrt(((3.5*((energy/100.)**(-0.8)))**2)+(0.15**2))
#print(s_p_e)
i = 0
while i<len(theta_array_bis):
	theta_array_bis_x = theta_array_bis[i]/s_p_e
	theta_array_bis_x_arr.append(theta_array_bis_x)
	i = i + 1


def king(x, f, sigma, gamma, sigma_t, gamma_t):
	y = f*((1./(2.*math.pi*(sigma**2)))*(1.-(1./gamma))*(1.+(x**2/(2.*(sigma**2)*gamma)))**(-gamma))+(1-f)*((1./(2.*math.pi*(sigma_t**2)))*(1.-(1./gamma_t))*(1.+(x**2/(2.*(sigma_t**2)*gamma_t)))**(-gamma_t))
	#y = (1./(2.*math.pi*(sigma**2)))*(1.-(1./gamma))*(1.+(x**2/(2.*(sigma**2)*gamma)))**(-gamma)
	return y

p, cov = scipy.optimize.curve_fit(king, theta_array_bis_x_arr, norm_N_bis/max(norm_N_bis))
#print(p)
#print(cov)
f = p[0]
sigma = p[1]
gamma = p[2]
sigma_t = p[3]
gamma_t = p[4]
y_arr = []
i = 0
while i<len(theta_array_bis_x_arr):
	y = king(theta_array_bis_x_arr[i],f, sigma, gamma, sigma_t, gamma_t)
	y_arr.append(y)
	i = i +1

y_arr = np.array(y_arr)

gdl = len(y_arr)-2
chi_quadro = (norm_N_bis[0]-y_arr[0])**2/y_arr[0]
i = 1 
while i < len(y_arr):
	chi_quadro = ((norm_N_bis[i]-y_arr[i])**2/y_arr[i]) + chi_quadro
	i = i + 1
chi_quadro_rid = chi_quadro/gdl
#print(chi_quadro_rid)

#chi_quadro = scipy.stats.chisquare(norm_N_bis, y_arr, len(y_arr)-2)

#print(chi_quadro)

err_y = err_norm_N_bis/max(err_norm_N_bis*2)
#err_y = (max(err_norm_N_bis)-err_norm_N_bis)/max(err_norm_N_bis)
fig = plt.figure(1, figsize=(10, 6))
ax = fig.add_subplot(111)

ax.errorbar(theta_array_bis, norm_N_bis/max(norm_N_bis), label = 'Dati simulati',  xerr= err_theta_bis, yerr = err_y, linestyle = 'None', color = 'red')  
#norm_N/norm_N[0] yerr = err_norm_N_bis
ax.plot(theta_array_bis, y_arr, label = 'Modello di King',color = 'blue', linewidth = 2.5, linestyle = 'dotted')
ax.set_yscale("log")	# set y log scale
#ax.set_xscale("log")	# set x log scale
ax.legend(loc = 'upper right', fontsize = 11)  
ax.set_ylim(0.0005,5.)
ax.set_title("eASTROGAM, E=100 MeV, $\\theta$ = 30$^\\circ$, $\\phi$ = 225$^\\circ$, $10^5$ fotoni")
ax.set_ylabel('Brillanza superficiale normalizzata [sr$^{-1}$]')    # x label
ax.set_xlabel('Distanza sferica [$^\\circ$]')	  # y label

plt.show()
