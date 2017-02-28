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

def strip_analysis(Total_vol_x_top, Total_vol_y_top, N_trig, Arch_vol_id_x_top, Arch_moth_id_x_top, Arch_Strip_id_x_top, Arch_Si_id_x_top, Arch_tray_id_x_top, Arch_plane_id_x_top, Arch_xpos_x_top, Arch_zpos_x_top, Arch_energy_dep_x_top, Arch_pair_flag_x_top, Arch_vol_id_y_top, Arch_moth_id_y_top, Arch_Strip_id_y_top, Arch_Si_id_y_top, Arch_tray_id_y_top, Arch_plane_id_y_top, Arch_ypos_y_top, Arch_zpos_y_top, Arch_energy_dep_y_top, Arch_pair_flag_y_top, event_id_tot, vol_id_tot, moth_id_tot, Strip_id_tot, Si_id_tot, tray_id_tot, plane_id_tot, energy_dep_tot, pair_flag_tot):

	Glob_event_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
	Glob_vol_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
	Glob_moth_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
	Glob_Strip_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
	Glob_Si_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
	Glob_tray_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
	Glob_plane_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
	Glob_xpos_x_top = np.zeros((Total_vol_x_top, N_trig))
	Glob_zpos_x_top = np.zeros((Total_vol_x_top, N_trig))
	Glob_energy_dep_x_top = np.zeros((Total_vol_x_top, N_trig))
	Glob_pair_flag_x_top = np.zeros((Total_vol_x_top, N_trig))

	Glob_event_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
	Glob_vol_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
	Glob_moth_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
	Glob_Strip_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
	Glob_Si_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
	Glob_tray_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
	Glob_plane_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
	Glob_ypos_y_top = np.zeros((Total_vol_y_top, N_trig))
	Glob_zpos_y_top = np.zeros((Total_vol_y_top, N_trig))
	Glob_energy_dep_y_top = np.zeros((Total_vol_y_top, N_trig))
	Glob_pair_flag_y_top = np.zeros((Total_vol_y_top, N_trig))



	for i in range(N_trig):
		Glob_vol_id_x_top[:,i] = Arch_vol_id_x_top
		Glob_moth_id_x_top[:,i] = Arch_moth_id_x_top
		Glob_Strip_id_x_top[:,i] = Arch_Strip_id_x_top
		Glob_Si_id_x_top[:,i] = Arch_Si_id_x_top
		Glob_tray_id_x_top[:,i] = Arch_tray_id_x_top
		Glob_plane_id_x_top[:,i] = Arch_plane_id_x_top
		Glob_xpos_x_top[:,i] = Arch_xpos_x_top
		Glob_zpos_x_top[:,i] = Arch_zpos_x_top
		Glob_energy_dep_x_top[:,i] = Arch_energy_dep_x_top
		Glob_pair_flag_x_top[:,i] = Arch_pair_flag_x_top

		Glob_vol_id_y_top[:,i] = Arch_vol_id_y_top
		Glob_moth_id_y_top[:,i] = Arch_moth_id_y_top
		Glob_Strip_id_y_top[:,i] = Arch_Strip_id_y_top
		Glob_Si_id_y_top[:,i] = Arch_Si_id_y_top
		Glob_tray_id_y_top[:,i] = Arch_tray_id_y_top
		Glob_plane_id_y_top[:,i] = Arch_plane_id_y_top
		Glob_ypos_y_top[:,i] = Arch_ypos_y_top
		Glob_zpos_y_top[:,i] = Arch_zpos_y_top
		Glob_energy_dep_y_top[:,i] = Arch_energy_dep_y_top
		Glob_pair_flag_y_top[:,i] = Arch_pair_flag_y_top


	j = 0
	N_ev = 0
	while 1:
		where_event_eq = np.where(event_id_tot == event_id_tot[j])
		where_event_eq = where_event_eq[0]

				
		event_id_temp = event_id_tot[where_event_eq]
		vol_id_temp = vol_id_tot[where_event_eq]
		moth_id_temp = moth_id_tot[where_event_eq]
		Strip_id_temp = Strip_id_tot[where_event_eq]
		Si_id_temp = Si_id_tot[where_event_eq]
		tray_id_temp = tray_id_tot[where_event_eq]
		plane_id_temp = plane_id_tot[where_event_eq]
		energy_dep_temp = energy_dep_tot[where_event_eq]
		pair_flag_temp = pair_flag_tot[where_event_eq]			


		# Funzione di ordinamento 1
										
		vol_sort_arr = np.argsort(vol_id_temp, kind='mergesort')
								
		vol_id_temp = vol_id_temp[vol_sort_arr]
		moth_id_temp = moth_id_temp[vol_sort_arr]
		Strip_id_temp = Strip_id_temp[vol_sort_arr]
		Si_id_temp = Si_id_temp[vol_sort_arr]
		tray_id_temp = tray_id_temp[vol_sort_arr]
		plane_id_temp = plane_id_temp[vol_sort_arr]
		energy_dep_temp = energy_dep_temp[vol_sort_arr]
		pair_flag_temp = pair_flag_temp[vol_sort_arr]				
				
			
		for z in range(Total_vol_x_top):
			where_hit_x_top_arr = np.where((Si_id_temp == 0) & (vol_id_temp == Glob_vol_id_x_top[z, N_ev]) & (moth_id_temp == Glob_moth_id_x_top[z, N_ev]))
			where_hit_x_top = where_hit_x_top_arr[0]

						
			if len(where_hit_x_top) != 0: 	
				Glob_energy_dep_x_top[z, N_ev] = energy_dep_temp[where_hit_x_top]
				Glob_pair_flag_x_top[z, N_ev] = pair_flag_temp[where_hit_x_top]
						
			where_hit_y_top_old = np.where((Si_id_temp == 1) & (vol_id_temp == Glob_vol_id_y_top[z, N_ev]) & (moth_id_temp == Glob_moth_id_y_top[z, N_ev]))
			where_hit_y_top = where_hit_y_top_old[0]

				
			if len(where_hit_y_top) != 0:
				Glob_energy_dep_y_top[z, N_ev] = energy_dep_temp[where_hit_y_top]
				Glob_pair_flag_y_top[z, N_ev] = pair_flag_temp[where_hit_y_top]

					

		N_event_eq = len(where_event_eq)
		if where_event_eq[N_event_eq-1] < len(event_id_tot)-1:
			j = where_event_eq[N_event_eq-1]+1
			N_ev = N_ev + 1
		else:
			break
			

	return(Glob_vol_id_x_top, Glob_moth_id_x_top, Glob_Strip_id_x_top, Glob_Si_id_x_top, Glob_tray_id_x_top, Glob_plane_id_x_top, Glob_xpos_x_top, Glob_zpos_x_top, Glob_energy_dep_x_top, Glob_pair_flag_x_top, Glob_vol_id_y_top, Glob_moth_id_y_top, Glob_Strip_id_y_top, Glob_Si_id_y_top, Glob_tray_id_y_top, Glob_plane_id_y_top, Glob_ypos_y_top, Glob_zpos_y_top, Glob_energy_dep_y_top, Glob_pair_flag_y_top, N_ev) 						

