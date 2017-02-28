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


def flag_events(event_id_tr, vol_id_tr, moth_id_tr, Strip_id_x, Strip_id_y, tray_id, plane_id, en_dep_tr, child_id_tr, proc_id_tr):

	e_dep_temp = [] 		
	event_id_tot = []
	vol_id_tot = []
	moth_id_tot = []
	Strip_id_x_tot = []
	Strip_id_y_tot = []
	tray_id_tot = []
	plane_id_tot = []
	energy_dep_tot = []
	pair_flag_tot = []
			
	j = 0
	while 1:          

		where_event_eq = np.where(event_id_tr == event_id_tr[j])
		where_event_eq = where_event_eq[0]
				
		vol_id_temp = vol_id_tr[where_event_eq]
		moth_id_temp  = moth_id_tr[where_event_eq]
		Strip_id_x_temp  = Strip_id_x[where_event_eq]
		Strip_id_y_temp  = Strip_id_y[where_event_eq]
		tray_id_temp  = tray_id[where_event_eq]
		plane_id_temp  = plane_id[where_event_eq]
		energy_dep_temp = en_dep_tr[where_event_eq]
		child_id_temp = child_id_tr[where_event_eq]
		proc_id_temp = proc_id_tr[where_event_eq]



		r = 0												
		while 1:
					
			where_vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]))										
			where_vol_eq = where_vol_eq[0]	

			where_other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]))
			where_other_vol = where_other_vol[0]

			e_dep_temp_old = np.sum(energy_dep_temp[where_vol_eq])
			event_id_tot_old = event_id_tr[j]
			vol_id_tot_old = vol_id_temp[r]
			moth_id_tot_old = moth_id_temp[r]
			Strip_id_x_tot_old = Strip_id_x_temp[r]
			Strip_id_y_tot_old = Strip_id_y_temp[r]
			tray_id_tot_old = tray_id_temp[r]
			plane_id_tot_old = plane_id_temp[r]
			energy_dep_tot_old = e_dep_temp_old


			e_dep_temp.append(e_dep_temp_old)
			event_id_tot.append(event_id_tot_old)
			vol_id_tot.append(vol_id_tot_old)
			moth_id_tot.append(moth_id_tot_old)
			Strip_id_x_tot.append(Strip_id_x_tot_old)
			Strip_id_y_tot.append(Strip_id_y_tot_old)
			tray_id_tot.append(tray_id_tot_old)
			plane_id_tot.append(plane_id_tot_old)
			energy_dep_tot.append(energy_dep_tot_old)					
					
			# if one of hits in the same volume is a pair the summed event is flagged as 1
			all_child = child_id_temp[where_vol_eq]
			all_proc = proc_id_temp[where_vol_eq]

			where_pair = np.where((all_child == 1) & (all_proc == 7))
			where_pair = where_pair[0]

			if len(where_pair) != 0:
				pair_flag_tot_old = 1
			else:
				pair_flag_tot_old = 0				
			pair_flag_tot.append(pair_flag_tot_old)	
 
			where_compton = np.where((all_child == 1) & (all_proc == 3))
			where_compton = where_compton[0]

			if len(where_compton) != 0:
				pair_flag_tot_old = 2
			else:
				pair_flag_tot_old = 0
				
			pair_flag_tot.append(pair_flag_tot_old)	

			if len(where_other_vol) != 0:
				vol_id_temp = vol_id_temp[where_other_vol]
				moth_id_temp = moth_id_temp[where_other_vol]
				Strip_id_x_temp = Strip_id_x_temp[where_other_vol]
				Strip_id_y_temp = Strip_id_y_temp[where_other_vol]
				tray_id_temp = tray_id_temp[where_other_vol]
				plane_id_temp = plane_id_temp[where_other_vol]
				energy_dep_temp = energy_dep_temp[where_other_vol]
				child_id_temp = child_id_temp[where_other_vol]
				proc_id_temp = proc_id_temp[where_other_vol]
					
			else:
				break
			
		N_event_eq = len(where_event_eq)                                          
		if where_event_eq[N_event_eq-1] < len(event_id_tr)-1:
			j = where_event_eq[N_event_eq-1]+1					
		else:
			break

	e_dep_temp = np.array(e_dep_temp) 		
	event_id_tot = np.array(event_id_tot)
	vol_id_tot = np.array(vol_id_tot)
	moth_id_tot = np.array(moth_id_tot)
	Strip_id_x_tot = np.array(Strip_id_x_tot)
	Strip_id_y_tot = np.array(Strip_id_y_tot)
	tray_id_tot = np.array(tray_id_tot)
	plane_id_tot = np.array(plane_id_tot)
	energy_dep_tot = np.array(energy_dep_tot)
	pair_flag_tot = np.array(pair_flag_tot)

	return(e_dep_temp, event_id_tot, vol_id_tot, moth_id_tot, Strip_id_x_tot, Strip_id_y_tot, tray_id_tot, plane_id_tot, energy_dep_tot, pair_flag_tot)
