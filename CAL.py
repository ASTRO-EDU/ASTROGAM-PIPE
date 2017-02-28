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


def G4_cal(event_id_cal, vol_id_cal, moth_id_cal, bar_ene, E_th_cal, cal_vol_start):

	N_trig_cal = 0

	event_id_tot_cal = []
	vol_id_tot_cal = []
	bar_id_tot = []
	moth_id_tot_cal = []
	bar_ene_tot = []
		
	j=0
	while 1:
			
		where_event_eq = np.where(event_id_cal == event_id_cal[j])
		where_event_eq = where_event_eq[0]

		N_trig_cal = N_trig_cal + 1
			
		vol_id_temp_cal = vol_id_cal[where_event_eq]
		moth_id_temp_cal = moth_id_cal[where_event_eq]
		bar_ene_temp = bar_ene[where_event_eq]

		r = 0
		while 1:
				 
			where_vol_eq = np.where(vol_id_temp_cal == vol_id_temp_cal[r])
			where_vol_eq = where_vol_eq[0]

			where_other_vol = np.where(vol_id_temp_cal != vol_id_temp_cal[r])
			where_other_vol = where_other_vol[0]

			bar_ene_tot_temp = np.sum(bar_ene_temp[where_vol_eq])
 
			if bar_ene_tot_temp >= E_th_cal:
				event_id_tot_cal_old = event_id_cal[j]
				vol_id_tot_cal_old = vol_id_temp_cal[r]
				bar_id_tot_old = vol_id_temp_cal[r] - cal_vol_start
				moth_id_tot_cal_old = moth_id_temp_cal[r]
				bar_ene_tot_old = np.sum(bar_ene_temp[where_vol_eq])

				event_id_tot_cal.append(event_id_tot_cal_old)
				vol_id_tot_cal.append(vol_id_tot_cal_old)
				bar_id_tot.append(bar_id_tot_old )
				moth_id_tot_cal.append(moth_id_tot_cal_old)
				bar_ene_tot.append(bar_ene_tot_old) 
					
			if where_other_vol != []:
				vol_id_temp_cal = vol_id_temp_cal[where_other_vol]
				moth_id_temp_cal = moth_id_temp_cal[where_other_vol]
				bar_ene_temp = bar_ene_temp[where_other_vol]
			else:
				break

			
		N_event_eq = len(where_event_eq)                                          
		if where_event_eq[N_event_eq-1] < len(event_id_cal)-1:
			j = where_event_eq[N_event_eq-1]+1
		else:
			break
		
		
	event_id_tot_cal = np.array((event_id_tot_cal), dtype = np.int64)
	bar_id_tot = np.array((bar_id_tot), dtype = np.int64)
	bar_ene_tot = np.array(bar_ene_tot)

	return(event_id_tot_cal, bar_id_tot, bar_ene_tot)
	
################################

def cal_sum(event_id_tot_cal, bar_ene_tot):

	event_id_tot_cal_sum = []
	bar_ene_tot_sum = []


	j=0
	while 1:
			
		where_event_eq = np.where(event_id_tot_cal == event_id_tot_cal[j])
		where_event_eq = where_event_eq[0]
			
			
		event_id_tot_cal_sum_old = event_id_tot_cal[j]
			
		event_id_tot_cal_sum.append(event_id_tot_cal_sum_old)
			
		bar_ene_tot_sum_old = np.sum(bar_ene_tot[where_event_eq])
		bar_ene_tot_sum = np.append(bar_ene_tot_sum, bar_ene_tot_sum_old)

		N_event_eq = len(where_event_eq)                                          
		if where_event_eq[N_event_eq-1] < len(event_id_tot_cal)-1:
			j = where_event_eq[N_event_eq-1]+1
		else:
			break


	event_id_tot_cal_sum = np.array((event_id_tot_cal_sum), dtype=np.int64)
	bar_ene_tot_sum = np.array(bar_ene_tot_sum)
	
	return(event_id_tot_cal_sum, bar_ene_tot_sum)	

