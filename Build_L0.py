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


def build_L0(Glob_energy_dep_x_top, Glob_energy_dep_y_top, Glob_vol_id_x_top, Glob_vol_id_y_top, Glob_moth_id_x_top, Glob_moth_id_y_top, Glob_Strip_id_x_top, Glob_Strip_id_y_top, Glob_Si_id_x_top, Glob_Si_id_y_top, Glob_tray_id_x_top, Glob_tray_id_y_top, Glob_plane_id_x_top, Glob_plane_id_y_top, Glob_xpos_x_top, Glob_ypos_y_top, Glob_zpos_x_top, Glob_zpos_y_top, Glob_pair_flag_x_top, Glob_pair_flag_y_top, N_trig, event_array):

	Glob_event_id_test = []
	Glob_vol_id_test = []
	Glob_moth_id_test = []
	Glob_Strip_id_test = []
	Glob_Si_id_test = []
	Glob_tray_id_test = []
	Glob_plane_id_test = []
	Glob_pos_test = []
	Glob_zpos_test = []
	Glob_energy_dep_test = []
	Glob_pair_flag_test = []

	for j in range(N_trig):

		where_test_x = np.where(Glob_energy_dep_x_top[:,j] > 0.)
		where_test_x = where_test_x[0]
				
		if len(where_test_x) != 0:
			Glob_vol_id_x_test_temp = Glob_vol_id_x_top[where_test_x,j]
			Glob_moth_id_x_test_temp = Glob_moth_id_x_top[where_test_x,j]
			Glob_Strip_id_x_test_temp = Glob_Strip_id_x_top[where_test_x,j]
			Glob_Si_id_x_test_temp = Glob_Si_id_x_top[where_test_x,j]
			Glob_tray_id_x_test_temp = Glob_tray_id_x_top[where_test_x,j]
			Glob_plane_id_x_test_temp = Glob_plane_id_x_top[where_test_x,j]
			Glob_xpos_x_test_temp = Glob_xpos_x_top[where_test_x,j]
			Glob_zpos_x_test_temp = Glob_zpos_x_top[where_test_x,j]
			Glob_energy_dep_x_test_temp = Glob_energy_dep_x_top[where_test_x,j]
			Glob_pair_flag_x_test_temp = Glob_pair_flag_x_top[where_test_x,j]

		where_test_y = np.where(Glob_energy_dep_y_top[:,j] > 0.)
		where_test_y = where_test_y[0]
				
		if len(where_test_y) != 0:
			Glob_vol_id_y_test_temp = Glob_vol_id_y_top[where_test_y,j]
			Glob_moth_id_y_test_temp = Glob_moth_id_y_top[where_test_y,j]
			Glob_Strip_id_y_test_temp = Glob_Strip_id_y_top[where_test_y,j]
			Glob_Si_id_y_test_temp = Glob_Si_id_y_top[where_test_y,j]
			Glob_tray_id_y_test_temp = Glob_tray_id_y_top[where_test_y,j]
			Glob_plane_id_y_test_temp = Glob_plane_id_y_top[where_test_y,j]
			Glob_ypos_y_test_temp = Glob_ypos_y_top[where_test_y,j]
			Glob_zpos_y_test_temp = Glob_zpos_y_top[where_test_y,j]
			Glob_energy_dep_y_test_temp = Glob_energy_dep_y_top[where_test_y,j]
			Glob_pair_flag_y_test_temp = Glob_pair_flag_y_top[where_test_y,j]
					
		if len(where_test_y) != 0 and len(where_test_x) != 0:
			Glob_vol_id_test_temp = np.concatenate((Glob_vol_id_y_test_temp, Glob_vol_id_x_test_temp))
			Glob_moth_id_test_temp = np.concatenate((Glob_moth_id_y_test_temp, Glob_moth_id_x_test_temp))
			Glob_Strip_id_test_temp = np.concatenate((Glob_Strip_id_y_test_temp, Glob_Strip_id_x_test_temp))
			Glob_Si_id_test_temp = np.concatenate((Glob_Si_id_y_test_temp, Glob_Si_id_x_test_temp))
			Glob_tray_id_test_temp = np.concatenate((Glob_tray_id_y_test_temp, Glob_tray_id_x_test_temp))
			Glob_plane_id_test_temp = np.concatenate((Glob_plane_id_y_test_temp, Glob_plane_id_x_test_temp))
			Glob_pos_test_temp = np.concatenate((Glob_ypos_y_test_temp, Glob_xpos_x_test_temp))
			Glob_zpos_test_temp = np.concatenate((Glob_zpos_y_test_temp, Glob_zpos_x_test_temp))
			Glob_energy_dep_test_temp = np.concatenate((Glob_energy_dep_y_test_temp, Glob_energy_dep_x_test_temp))
			Glob_pair_flag_test_temp = np.concatenate((Glob_pair_flag_y_test_temp, Glob_pair_flag_x_test_temp))

		elif len(where_test_y) != 0 and len(where_test_x) == 0:
			Glob_vol_id_test_temp = Glob_vol_id_y_test_temp
			Glob_moth_id_test_temp = Glob_moth_id_y_test_temp
			Glob_Strip_id_test_temp = Glob_Strip_id_y_test_temp
			Glob_Si_id_test_temp = Glob_Si_id_y_test_temp
			Glob_tray_id_test_temp = Glob_tray_id_y_test_temp
			Glob_plane_id_test_temp = Glob_plane_id_y_test_temp
			Glob_pos_test_temp = Glob_ypos_y_test_temp
			Glob_zpos_test_temp = Glob_zpos_y_test_temp
			Glob_energy_dep_test_temp = Glob_energy_dep_y_test_temp
			Glob_pair_flag_test_temp = Glob_pair_flag_y_test_temp
          				
		elif len(where_test_y) == 0 and len(where_test_x) != 0:
			Glob_vol_id_test_temp = Glob_vol_id_x_test_temp
			Glob_moth_id_test_temp = Glob_moth_id_x_test_temp
			Glob_Strip_id_test_temp = Glob_Strip_id_x_test_temp
			Glob_Si_id_test_temp = Glob_Si_id_x_test_temp
			Glob_tray_id_test_temp = Glob_tray_id_x_test_temp
			Glob_plane_id_test_temp = Glob_plane_id_x_test_temp
			Glob_pos_test_temp = Glob_xpos_x_test_temp
			Glob_zpos_test_temp = Glob_zpos_x_test_temp
			Glob_energy_dep_test_temp = Glob_energy_dep_x_test_temp
			Glob_pair_flag_test_temp = Glob_pair_flag_x_test_temp

		# Funzione di ordinamento 2
										
		tray_sort_arr_temp = np.argsort(Glob_tray_id_test_temp, kind='mergesort')				

		tray_sort_arr = tray_sort_arr_temp[::-1]
		Glob_vol_id_test_temp = Glob_vol_id_test_temp[tray_sort_arr]
	        Glob_moth_id_test_temp = Glob_moth_id_test_temp[tray_sort_arr]
		Glob_Strip_id_test_temp = Glob_Strip_id_test_temp[tray_sort_arr]
		Glob_Si_id_test_temp = Glob_Si_id_test_temp[tray_sort_arr]
		Glob_tray_id_test_temp = Glob_tray_id_test_temp[tray_sort_arr]
		Glob_plane_id_test_temp = Glob_plane_id_test_temp[tray_sort_arr]
		Glob_pos_test_temp = Glob_pos_test_temp[tray_sort_arr]
		Glob_zpos_test_temp = Glob_zpos_test_temp[tray_sort_arr]
		Glob_energy_dep_test_temp = Glob_energy_dep_test_temp[tray_sort_arr]
		Glob_pair_flag_test_temp = Glob_pair_flag_test_temp[tray_sort_arr]


		vol_id_intray = []
		moth_id_intray = []
		Strip_id_intray = []
		Si_id_intray = []
		tray_id_intray = []
		plane_id_intray = []
		pos_intray = []
		zpos_intray = []
		energy_dep_intray = []
		pair_flag_intray = []

		
		intray = 0
		while 1:
			where_tray_eq = np.where(Glob_tray_id_test_temp == Glob_tray_id_test_temp[intray])
			where_tray_eq = where_tray_eq[0]

			where_other_tray = np.where(Glob_tray_id_test_temp != Glob_tray_id_test_temp[intray])
			where_other_tray = where_other_tray[0]


			vol_id_extract = Glob_vol_id_test_temp[where_tray_eq]
			moth_id_extract = Glob_moth_id_test_temp[where_tray_eq]
			Strip_id_extract = Glob_Strip_id_test_temp[where_tray_eq]
			Si_id_extract = Glob_Si_id_test_temp[where_tray_eq]
			tray_id_extract = Glob_tray_id_test_temp[where_tray_eq]
			plane_id_extract = Glob_plane_id_test_temp[where_tray_eq]
			pos_extract = Glob_pos_test_temp[where_tray_eq]
			zpos_extract = Glob_zpos_test_temp[where_tray_eq]
			energy_dep_extract = Glob_energy_dep_test_temp[where_tray_eq]
			pair_flag_extract = Glob_pair_flag_test_temp[where_tray_eq]

			where_Y = np.where(Si_id_extract == 1)
			where_Y = where_Y[0] 

			if len(where_Y) != 0:
				vol_id_intray_old = vol_id_extract[where_Y]
				moth_id_intray_old = moth_id_extract[where_Y]
				Strip_id_intray_old = Strip_id_extract[where_Y]
				Si_id_intray_old = Si_id_extract[where_Y]
				tray_id_intray_old = tray_id_extract[where_Y]
				plane_id_intray_old = plane_id_extract[where_Y]
				pos_intray_old = pos_extract[where_Y]
				zpos_intray_old = zpos_extract[where_Y]
				energy_dep_intray_old = energy_dep_extract[where_Y]
				pair_flag_intray_old = pair_flag_extract[where_Y]

				vol_id_intray.append(vol_id_intray_old)
				moth_id_intray.append(moth_id_intray_old)
				Strip_id_intray.append(Strip_id_intray_old)
				Si_id_intray.append(Si_id_intray_old)
				tray_id_intray.append(tray_id_intray_old)
				plane_id_intray.append(plane_id_intray_old)
				pos_intray.append(pos_intray_old)
				zpos_intray.append(zpos_intray_old)
				energy_dep_intray.append(energy_dep_intray_old)
				pair_flag_intray.append(pair_flag_intray_old)


			where_X = np.where(Si_id_extract == 0)
			where_X = where_X[0]

			if len(where_X) != 0:

				vol_id_intray_old = vol_id_extract[where_X]
				moth_id_intray_old = moth_id_extract[where_X]
				Strip_id_intray_old = Strip_id_extract[where_X]
				Si_id_intray_old = Si_id_extract[where_X]
				tray_id_intray_old = tray_id_extract[where_X]
				plane_id_intray_old = plane_id_extract[where_X]
				pos_intray_old = pos_extract[where_X]
				zpos_intray_old = zpos_extract[where_X]
				energy_dep_intray_old = energy_dep_extract[where_X]
				pair_flag_intray_old = pair_flag_extract[where_X]

				vol_id_intray.append(vol_id_intray_old)
				moth_id_intray.append(moth_id_intray_old)
				Strip_id_intray.append(Strip_id_intray_old)
				Si_id_intray.append(Si_id_intray_old)
				tray_id_intray.append(tray_id_intray_old)
				plane_id_intray.append(plane_id_intray_old)
				pos_intray.append(pos_intray_old)
				zpos_intray.append(zpos_intray_old)
				energy_dep_intray.append(energy_dep_intray_old)
				pair_flag_intray.append(pair_flag_intray_old)
 
			N_tray_eq = len(where_tray_eq)
			if where_tray_eq[N_tray_eq-1] < len(Glob_tray_id_test_temp)-1:
				intray = where_tray_eq[N_tray_eq-1]+1						
			else:
				break 
 
		vol_id_intray = np.ma.concatenate(vol_id_intray)
		moth_id_intray = np.ma.concatenate(moth_id_intray)
		Strip_id_intray = np.ma.concatenate(Strip_id_intray)
		Si_id_intray = np.ma.concatenate(Si_id_intray)
		tray_id_intray = np.ma.concatenate(tray_id_intray)
		plane_id_intray = np.ma.concatenate(plane_id_intray)
		pos_intray = np.ma.concatenate(pos_intray)
		zpos_intray = np.ma.concatenate(zpos_intray)
		energy_dep_intray = np.ma.concatenate(energy_dep_intray)
		pair_flag_intray = np.ma.concatenate(pair_flag_intray)
 				

		vol_id_temp = np.array(vol_id_intray)
		moth_id_temp = np.array(moth_id_intray)
		Strip_id_temp = np.array(Strip_id_intray)
		Si_id_temp = np.array(Si_id_intray)
		tray_id_temp = np.array(tray_id_intray)
		plane_id_temp = np.array(plane_id_intray)
		pos_temp = np.array(pos_intray)
		zpos_temp = np.array(zpos_intray)
		energy_dep_temp = np.array(energy_dep_intray)
		pair_flag_temp = np.array(pair_flag_intray)				  
 
 
		event_id_temp = []
		for k in range(len(vol_id_temp)):
			event_id_temp_old = event_array[j]
			event_id_temp.append(event_id_temp_old)


		event_id_temp = np.array(event_id_temp)
				
				
		Glob_event_id_test.append(event_id_temp)
		Glob_vol_id_test.append(vol_id_temp)
		Glob_moth_id_test.append(moth_id_temp)
		Glob_Strip_id_test.append(Strip_id_temp)
		Glob_Si_id_test.append(Si_id_temp)
		Glob_tray_id_test.append(tray_id_temp)
		Glob_plane_id_test.append(plane_id_temp)
		Glob_pos_test.append(pos_temp)
		Glob_zpos_test.append(zpos_temp)
		Glob_energy_dep_test.append(energy_dep_temp)
		Glob_pair_flag_test.append(pair_flag_temp)


	Glob_event_id_test = np.ma.concatenate(Glob_event_id_test)	
	Glob_vol_id_test = np.ma.concatenate(Glob_vol_id_test)
	Glob_moth_id_test = np.ma.concatenate(Glob_moth_id_test)
	Glob_Strip_id_test = np.ma.concatenate(Glob_Strip_id_test)
	Glob_Si_id_test = np.ma.concatenate(Glob_Si_id_test)
	Glob_tray_id_test = np.ma.concatenate(Glob_tray_id_test)
	Glob_plane_id_test = np.ma.concatenate(Glob_plane_id_test)
	Glob_pos_test = np.ma.concatenate(Glob_pos_test)
	Glob_zpos_test = np.ma.concatenate(Glob_zpos_test)
	Glob_energy_dep_test = np.ma.concatenate(Glob_energy_dep_test)
	Glob_pair_flag_test = np.ma.concatenate(Glob_pair_flag_test)			

	Glob_event_id_test = np.array((Glob_event_id_test), dtype=np.int64)
	Glob_vol_id_test = np.array((Glob_vol_id_test), dtype=np.int64)
	Glob_moth_id_test = np.array((Glob_moth_id_test), dtype=np.int64)
	Glob_Strip_id_test = np.array((Glob_Strip_id_test), dtype=np.int64)
	Glob_Si_id_test = np.array((Glob_Si_id_test), dtype=np.int64)
	Glob_tray_id_test = np.array((Glob_tray_id_test), dtype=np.int64)
	Glob_plane_id_test = np.array((Glob_plane_id_test), dtype=np.int64)
	Glob_pos_test = np.array(Glob_pos_test)
	Glob_zpos_test = np.array(Glob_zpos_test)
	Glob_energy_dep_test = np.array(Glob_energy_dep_test)
	Glob_pair_flag_test = np.array((Glob_pair_flag_test), dtype=np.int64)

	return(Glob_event_id_test, Glob_vol_id_test, Glob_moth_id_test, Glob_Strip_id_test, Glob_Si_id_test, Glob_tray_id_test, Glob_plane_id_test, Glob_pos_test, Glob_zpos_test, Glob_energy_dep_test, Glob_pair_flag_test)    
