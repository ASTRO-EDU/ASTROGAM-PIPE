from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import pickle
from astropy.table import Table
import copy
import time
import sys

# Summing the energy along the strip

def summing_energy(E_th,event_id_tot_temp, vol_id_tot_temp, moth_id_tot_temp, Strip_id_tot_temp, Si_id_tot_temp, tray_id_tot_temp, plane_id_tot_temp, energy_dep_tot_temp, pair_flag_tot_temp):

	e_dep_temp = [] 		
	event_id_tot = []
	vol_id_tot = []
	moth_id_tot = []
	Strip_id_tot = []
	Si_id_tot = []
	tray_id_tot = []
	plane_id_tot = []
	energy_dep_tot = []
	pair_flag_tot = []


	j = 0
	while 1:           
		where_event_eq = np.where(event_id_tot_temp == event_id_tot_temp[j])
		where_event_eq = where_event_eq[0]


		vol_id_temp = vol_id_tot_temp[where_event_eq]
		moth_id_temp = moth_id_tot_temp[where_event_eq]
		Strip_id_temp = Strip_id_tot_temp[where_event_eq]
		Si_id_temp = Si_id_tot_temp[where_event_eq]
		tray_id_temp = tray_id_tot_temp[where_event_eq]
		plane_id_temp = plane_id_tot_temp[where_event_eq]
		energy_dep_temp = energy_dep_tot_temp[where_event_eq]
		pair_flag_temp = pair_flag_tot_temp[where_event_eq]
				

				
		r = 0												
		while 1:

			where_vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]) & (Si_id_temp == 0))	
			where_vol_eq = where_vol_eq[0]
					
			where_other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]) | (Si_id_temp != 0))
			where_other_vol = where_other_vol[0]

			if len(where_vol_eq) != 0:
				e_dep_temp_old = np.sum(energy_dep_temp[where_vol_eq])
				event_id_tot_old = event_id_tot_temp[j]
				vol_id_tot_old = vol_id_temp[r]
				moth_id_tot_old = moth_id_temp[r]
				Si_id_tot_old = 0
				Strip_id_tot_old = Strip_id_temp[r]
				tray_id_tot_old = tray_id_temp[r]
				plane_id_tot_old = plane_id_temp[r]
				energy_dep_tot_old = e_dep_temp_old
				pair_flag_tot_old = pair_flag_temp[r]

				e_dep_temp.append(e_dep_temp_old)
				event_id_tot.append(event_id_tot_old)
				vol_id_tot.append(vol_id_tot_old)
				moth_id_tot.append(moth_id_tot_old)
				Strip_id_tot.append(Strip_id_tot_old)
				Si_id_tot.append(Si_id_tot_old)
				tray_id_tot.append(tray_id_tot_old)
				plane_id_tot.append(plane_id_tot_old)
				energy_dep_tot.append(energy_dep_tot_old)					
				pair_flag_tot.append(pair_flag_tot_old)

				if len(where_other_vol) != 0:
					vol_id_temp = vol_id_temp[where_other_vol]
					moth_id_temp = moth_id_temp[where_other_vol]
					Strip_id_temp = Strip_id_temp[where_other_vol]
					Si_id_temp = Si_id_temp[where_other_vol]
					tray_id_temp = tray_id_temp[where_other_vol]
					plane_id_temp = plane_id_temp[where_other_vol]
					energy_dep_temp = energy_dep_temp[where_other_vol]
					pair_flag_temp = pair_flag_temp[where_other_vol]
				else:							
					break


			where_vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]) & (Si_id_temp == 1))	
			where_vol_eq = where_vol_eq[0]
					
			where_other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]) | (Si_id_temp != 1))
			where_other_vol = where_other_vol[0]
					
			if len(where_vol_eq) != 0:
				e_dep_temp_old = np.sum(energy_dep_temp[where_vol_eq])
				event_id_tot_old = event_id_tot_temp[j]
				vol_id_tot_old = vol_id_temp[r]
				moth_id_tot_old = moth_id_temp[r]
				Si_id_tot_old = 1
				Strip_id_tot_old = Strip_id_temp[r]
				tray_id_tot_old = tray_id_temp[r]
				plane_id_tot_old = plane_id_temp[r]
				energy_dep_tot_old = e_dep_temp_old
				pair_flag_tot_old = pair_flag_temp[r]

				e_dep_temp.append(e_dep_temp_old)
				event_id_tot.append(event_id_tot_old)
				vol_id_tot.append(vol_id_tot_old)
				moth_id_tot.append(moth_id_tot_old)
				Strip_id_tot.append(Strip_id_tot_old)
				Si_id_tot.append(Si_id_tot_old)
				tray_id_tot.append(tray_id_tot_old)
				plane_id_tot.append(plane_id_tot_old)
				energy_dep_tot.append(energy_dep_tot_old)					
				pair_flag_tot.append(pair_flag_tot_old)

				if len(where_other_vol) != 0:
					vol_id_temp = vol_id_temp[where_other_vol]
					moth_id_temp = moth_id_temp[where_other_vol]
					Strip_id_temp = Strip_id_temp[where_other_vol]
					Si_id_temp = Si_id_temp[where_other_vol]
					tray_id_temp = tray_id_temp[where_other_vol]
					plane_id_temp = plane_id_temp[where_other_vol]
					energy_dep_temp = energy_dep_temp[where_other_vol]
					pair_flag_temp = pair_flag_temp[where_other_vol]
				else:
					break
			
		N_event_eq = len(where_event_eq)                                          
		if where_event_eq[N_event_eq-1] < len(event_id_tot_temp)-1:
			j = where_event_eq[N_event_eq-1]+1					
		else:
			break


	e_dep_temp = np.array(e_dep_temp) 		
	event_id_tot = np.array(event_id_tot)
	vol_id_tot = np.array(vol_id_tot)
	moth_id_tot = np.array(moth_id_tot)
	Strip_id_tot = np.array(Strip_id_tot)
	Si_id_tot = np.array(Si_id_tot)
	tray_id_tot = np.array(tray_id_tot)
	plane_id_tot = np.array(plane_id_tot)
	energy_dep_tot = np.array(energy_dep_tot)
	pair_flag_tot = np.array(pair_flag_tot)

	# apply the energy thresold
			
	where_eth = np.where(energy_dep_tot >= E_th)
	where_eth = where_eth[0]
			

	event_id_tot = event_id_tot[where_eth]
	vol_id_tot = vol_id_tot[where_eth]
	moth_id_tot = moth_id_tot[where_eth]
	Strip_id_tot = Strip_id_tot[where_eth]
	Si_id_tot = Si_id_tot[where_eth]
	tray_id_tot = tray_id_tot[where_eth]
	plane_id_tot = plane_id_tot[where_eth]
	energy_dep_tot = energy_dep_tot[where_eth]
	pair_flag_tot = pair_flag_tot[where_eth]


	return(event_id_tot, vol_id_tot, moth_id_tot, Strip_id_tot, Si_id_tot, tray_id_tot, plane_id_tot, energy_dep_tot, pair_flag_tot)  
