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


def AC_analysis(event_id_ac, vol_id_ac, energy_dep_ac, panel_S, panel_D, panel_F, panel_B, panel_top):

	N_trig_ac = 0


	event_id_tot_ac = []
	vol_id_tot_ac = []
	moth_id_tot_ac = []
	energy_dep_tot_ac = []



	j=0
	while 1:
		where_event_eq = np.where(event_id_ac == event_id_ac[j])
		where_event_eq = where_event_eq[0]

		N_trig_ac = N_trig_ac + 1

		vol_id_temp_ac = vol_id_ac[where_event_eq]
		moth_id_temp_ac = moth_id_ac[where_event_eq]
		energy_dep_temp_ac = energy_dep_ac[where_event_eq]


		r = 0
		while 1:
			where_vol_eq = np.where(vol_id_temp_ac == vol_id_temp_ac[r])
			where_vol_eq = where_vol_eq[0]

			where_other_vol = np.where(vol_id_temp_ac != vol_id_temp_ac[r])
			where_other_vol = where_other_vol[0]

			event_id_tot_ac_old = event_id_ac[j]
			vol_id_tot_ac_old = vol_id_temp_ac[r]
			moth_id_tot_ac_old = moth_id_temp_ac[r]
			energy_dep_tot_ac_old = np.sum(energy_dep_temp_ac[where_vol_eq])
	
			event_id_tot_ac.append(event_id_tot_ac_old)
			vol_id_tot_ac.append(vol_id_tot_ac_old)
			moth_id_tot_ac.append(moth_id_tot_ac_old)
			energy_dep_tot_ac.append(energy_dep_tot_ac_old)


			if len(where_other_vol) != 0:
				vol_id_temp_ac = vol_id_temp_ac(where_other_vol)
				moth_id_temp_ac = moth_id_temp_ac(where_other_vol)
				energy_dep_temp_ac = energy_dep_temp_ac(where_other_vol)
			else:
				break


		N_event_eq = len(where_event_eq)                                          
		if where_event_eq[N_event_eq-1] < len(event_id_ac)-1:
			j = where_event_eq[N_event_eq-1]+1
		else:
			break



	# AC panel IDs

	AC_panel = []
	AC_subpanel = []

	for j in range(len(vol_id_tot_ac)):

		if vol_id_tot_ac[j] >= panel_S[0] and vol_id_tot_ac[j] <= panel_S[2]:
			AC_panel_old[j] = 'S'
			AC_panel.append(AC_panel_old)

			if vol_id_tot_ac[j] == panel_S[0]:
				AC_subpanel_old[j] = 3
				AC_subpanel.append(AC_subpanel_old)
				
			if vol_id_tot_ac[j] == panel_S[1]:
				AC_subpanel_old[j] = 2
				AC_subpanel.append(AC_subpanel_old)

			if vol_id_tot_ac[j] == panel_S[2]:
				AC_subpanel_old[j] = 1
				AC_subpanel.append(AC_subpanel_old)


		if vol_id_tot_ac[j] >= panel_D[0] and vol_id_tot_ac[j] <= panel_D[2]:
			AC_panel_old[j] = 'D'
			AC_panel.append(AC_panel_old)

			if vol_id_tot_ac[j] == panel_D[0]:
				AC_subpanel_old[j] = 3
				AC_subpanel.append(AC_subpanel_old)
				
			if vol_id_tot_ac[j] == panel_D[1]:
				AC_subpanel_old[j] = 2
				AC_subpanel.append(AC_subpanel_old)

			if vol_id_tot_ac[j] == panel_D[2]:
				AC_subpanel_old[j] = 1
				AC_subpanel.append(AC_subpanel_old)


		if vol_id_tot_ac[j] >= panel_F[0] and vol_id_tot_ac[j] <= panel_F[2]:
			AC_panel_old[j] = 'F'
			AC_panel.append(AC_panel_old)

			if vol_id_tot_ac[j] == panel_F[0]:
				AC_subpanel_old[j] = 3
				AC_subpanel.append(AC_subpanel_old)
				
			if vol_id_tot_ac[j] == panel_F[1]:
				AC_subpanel_old[j] = 2
				AC_subpanel.append(AC_subpanel_old)

			if vol_id_tot_ac[j] == panel_F[2]:
				AC_subpanel_old[j] = 1
				AC_subpanel.append(AC_subpanel_old)


		if vol_id_tot_ac[j] >= panel_B[0] and vol_id_tot_ac[j] <= panel_B[2]:
			AC_panel_old[j] = 'B'
			AC_panel.append(AC_panel_old)

			if vol_id_tot_ac[j] == panel_B[0]:
				AC_subpanel_old[j] = 3
				AC_subpanel.append(AC_subpanel_old)
				
			if vol_id_tot_ac[j] == panel_B[1]:
				AC_subpanel_old[j] = 2
				AC_subpanel.append(AC_subpanel_old)

			if vol_id_tot_ac[j] == panel_B[2]:
				AC_subpanel_old[j] = 1
				AC_subpanel.append(AC_subpanel_old)


		if vol_id_tot_ac[j] == panel_top:
			AC_panel_old[j] = 'T'
			AC_subpanel_old[j] = 0

			AC_panel.append(AC_panel_old)
			AC_subpanel.append(AC_subpanel_old)


	
	AC_panel = np.array(AC_panel)
	AC_subpanel = np.array(AC_subpanel)
	event_id_tot_ac = np.array(event_id_tot_ac)
	energy_dep_tot_ac = np.array(energy_dep_tot_ac)

	return(AC_panel, AC_subpanel, event_id_tot_ac, energy_dep_tot_ac)

