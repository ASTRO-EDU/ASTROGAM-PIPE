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

def merging(event_array, N_trig, Glob_event_id_x_top_cluster, Glob_event_id_y_top_cluster, Glob_Strip_number_x_top_cluster, Glob_Strip_number_y_top_cluster, Glob_Si_id_x_top_cluster, Glob_Si_id_y_top_cluster, Glob_tray_id_x_top_cluster, Glob_tray_id_y_top_cluster, Glob_plane_id_x_top_cluster, Glob_plane_id_y_top_cluster, Glob_xpos_x_top_cluster, Glob_ypos_y_top_cluster, Glob_zpos_x_top_cluster, Glob_zpos_y_top_cluster, Glob_energy_dep_x_top_cluster, Glob_energy_dep_y_top_cluster, Glob_pair_flag_x_top_cluster, Glob_pair_flag_y_top_cluster):

	Glob_event_id_cluster = []
	Glob_Si_id_cluster = []
	Glob_tray_id_cluster = []
	Glob_plane_id_cluster = []
	Glob_pos_cluster = []
	Glob_zpos_cluster = []
	Glob_energy_dep_cluster = []
	Glob_Strip_number_cluster = []
	Glob_pair_flag_cluster = []
			
			
	for j in range(N_trig):
		Glob_Strip_number_cluster_temp = []
		Glob_Si_id_cluster_temp = []
		Glob_tray_id_cluster_temp = []
		Glob_plane_id_cluster_temp = []
		Glob_pos_cluster_temp = []
		Glob_zpos_cluster_temp = []
		Glob_energy_dep_cluster_temp = []
		Glob_pair_flag_cluster_temp = []

		where_cluster_x_top = np.where(Glob_event_id_x_top_cluster == j)
		where_cluster_x_top = where_cluster_x_top[0]

		if len(where_cluster_x_top) != 0:
					
			Glob_Strip_number_cluster_temp_old = Glob_Strip_number_x_top_cluster[where_cluster_x_top]
			Glob_Si_id_cluster_temp_old = Glob_Si_id_x_top_cluster[where_cluster_x_top]
			Glob_tray_id_cluster_temp_old = Glob_tray_id_x_top_cluster[where_cluster_x_top]
			Glob_plane_id_cluster_temp_old = Glob_plane_id_x_top_cluster[where_cluster_x_top]
			Glob_pos_cluster_temp_old = Glob_xpos_x_top_cluster[where_cluster_x_top]
			Glob_zpos_cluster_temp_old = Glob_zpos_x_top_cluster[where_cluster_x_top]
			Glob_energy_dep_cluster_temp_old = Glob_energy_dep_x_top_cluster[where_cluster_x_top]
			Glob_pair_flag_cluster_temp_old = Glob_pair_flag_x_top_cluster[where_cluster_x_top]

			Glob_Strip_number_cluster_temp = np.append(Glob_Strip_number_cluster_temp, Glob_Strip_number_cluster_temp_old)
			Glob_Si_id_cluster_temp = np.append(Glob_Si_id_cluster_temp, Glob_Si_id_cluster_temp_old)
			Glob_tray_id_cluster_temp = np.append(Glob_tray_id_cluster_temp, Glob_tray_id_cluster_temp_old)
			Glob_plane_id_cluster_temp = np.append(Glob_plane_id_cluster_temp, Glob_plane_id_cluster_temp_old)
			Glob_pos_cluster_temp = np.append(Glob_pos_cluster_temp, Glob_pos_cluster_temp_old)
			Glob_zpos_cluster_temp = np.append(Glob_zpos_cluster_temp, Glob_zpos_cluster_temp_old)
			Glob_energy_dep_cluster_temp = np.append(Glob_energy_dep_cluster_temp, Glob_energy_dep_cluster_temp_old)
			Glob_pair_flag_cluster_temp = np.append(Glob_pair_flag_cluster_temp, Glob_pair_flag_cluster_temp_old)


			where_cluster_y_top = np.where(Glob_event_id_y_top_cluster == j)
			where_cluster_y_top = where_cluster_y_top[0]

			if len(where_cluster_y_top) != 0:

				Glob_Strip_number_cluster_temp_old = Glob_Strip_number_y_top_cluster[where_cluster_y_top]
				Glob_Si_id_cluster_temp_old = Glob_Si_id_y_top_cluster[where_cluster_y_top]
				Glob_tray_id_cluster_temp_old = Glob_tray_id_y_top_cluster[where_cluster_y_top]
				Glob_plane_id_cluster_temp_old = Glob_plane_id_y_top_cluster[where_cluster_y_top]
				Glob_pos_cluster_temp_old = Glob_ypos_y_top_cluster[where_cluster_y_top]
				Glob_zpos_cluster_temp_old = Glob_zpos_y_top_cluster[where_cluster_y_top]
				Glob_energy_dep_cluster_temp_old = Glob_energy_dep_y_top_cluster[where_cluster_y_top]
				Glob_pair_flag_cluster_temp_old = Glob_pair_flag_y_top_cluster[where_cluster_y_top]

				Glob_Strip_number_cluster_temp = np.append(Glob_Strip_number_cluster_temp, Glob_Strip_number_cluster_temp_old)
				Glob_Si_id_cluster_temp = np.append(Glob_Si_id_cluster_temp, Glob_Si_id_cluster_temp_old)
				Glob_tray_id_cluster_temp = np.append(Glob_tray_id_cluster_temp, Glob_tray_id_cluster_temp_old)
				Glob_plane_id_cluster_temp = np.append(Glob_plane_id_cluster_temp, Glob_plane_id_cluster_temp_old)
				Glob_pos_cluster_temp = np.append(Glob_pos_cluster_temp, Glob_pos_cluster_temp_old)
				Glob_zpos_cluster_temp = np.append(Glob_zpos_cluster_temp, Glob_zpos_cluster_temp_old)
				Glob_energy_dep_cluster_temp = np.append(Glob_energy_dep_cluster_temp, Glob_energy_dep_cluster_temp_old)
				Glob_pair_flag_cluster_temp = np.append(Glob_pair_flag_cluster_temp, Glob_pair_flag_cluster_temp_old)

		else:
					
			where_cluster_y_top = np.where(Glob_event_id_y_top_cluster == j)
			where_cluster_y_top = where_cluster_y_top[0]

			if len(where_cluster_y_top) != 0:

				Glob_Strip_number_cluster_temp_old = Glob_Strip_number_y_top_cluster[where_cluster_y_top]
				Glob_Si_id_cluster_temp_old = Glob_Si_id_y_top_cluster[where_cluster_y_top]
				Glob_tray_id_cluster_temp_old = Glob_tray_id_y_top_cluster[where_cluster_y_top]
				Glob_plane_id_cluster_temp_old = Glob_plane_id_y_top_cluster[where_cluster_y_top]
				Glob_pos_cluster_temp_old = Glob_ypos_y_top_cluster[where_cluster_y_top]
				Glob_zpos_cluster_temp_old = Glob_zpos_y_top_cluster[where_cluster_y_top]
				Glob_energy_dep_cluster_temp_old = Glob_energy_dep_y_top_cluster[where_cluster_y_top]
				Glob_pair_flag_cluster_temp_old = Glob_pair_flag_y_top_cluster[where_cluster_y_top]

				Glob_Strip_number_cluster_temp = np.append(Glob_Strip_number_cluster_temp, Glob_Strip_number_cluster_temp_old)
				Glob_Si_id_cluster_temp = np.append(Glob_Si_id_cluster_temp, Glob_Si_id_cluster_temp_old)
				Glob_tray_id_cluster_temp = np.append(Glob_tray_id_cluster_temp, Glob_tray_id_cluster_temp_old)
				Glob_plane_id_cluster_temp = np.append(Glob_plane_id_cluster_temp, Glob_plane_id_cluster_temp_old)
				Glob_pos_cluster_temp = np.append(Glob_pos_cluster_temp, Glob_pos_cluster_temp_old)
				Glob_zpos_cluster_temp = np.append(Glob_zpos_cluster_temp, Glob_zpos_cluster_temp_old)
				Glob_energy_dep_cluster_temp = np.append(Glob_energy_dep_cluster_temp, Glob_energy_dep_cluster_temp_old)
				Glob_pair_flag_cluster_temp = np.append(Glob_pair_flag_cluster_temp, Glob_pair_flag_cluster_temp_old)
				
				
		# Funzione di ordinamento 7

		tray_sort_arr_temp = np.argsort(Glob_tray_id_cluster_temp, kind='mergesort')				

		# fine funzione di ordinamento	


		tray_sort_arr = tray_sort_arr_temp[::-1]
				
		Glob_Si_id_cluster_temp = Glob_Si_id_cluster_temp[tray_sort_arr]
		Glob_tray_id_cluster_temp = Glob_tray_id_cluster_temp[tray_sort_arr]
		Glob_plane_id_cluster_temp = Glob_plane_id_cluster_temp[tray_sort_arr]
		Glob_pos_cluster_temp = Glob_pos_cluster_temp[tray_sort_arr]
		Glob_zpos_cluster_temp = Glob_zpos_cluster_temp[tray_sort_arr]
		Glob_energy_dep_cluster_temp = Glob_energy_dep_cluster_temp[tray_sort_arr]
		Glob_Strip_number_cluster_temp = Glob_Strip_number_cluster_temp[tray_sort_arr]
		Glob_pair_flag_cluster_temp = Glob_pair_flag_cluster_temp[tray_sort_arr]
				
		Si_id_temp = []
		tray_id_temp = []
		plane_id_temp = []
		pos_temp = []
		zpos_temp = []
		energy_dep_temp = []
		strip_number_temp = []
		pair_flag_temp = []				
				
				
		intray = 0
		while 1:
			where_tray_eq = np.where(Glob_tray_id_cluster_temp == Glob_tray_id_cluster_temp[intray])
			where_tray_eq = where_tray_eq[0]
					
			Si_id_extract = Glob_Si_id_cluster_temp[where_tray_eq]
			tray_id_extract = Glob_tray_id_cluster_temp[where_tray_eq]
			plane_id_extract = Glob_plane_id_cluster_temp[where_tray_eq]
			pos_extract = Glob_pos_cluster_temp[where_tray_eq]
			zpos_extract = Glob_zpos_cluster_temp[where_tray_eq]
			energy_dep_extract = Glob_energy_dep_cluster_temp[where_tray_eq]
			strip_number_extract = Glob_Strip_number_cluster_temp[where_tray_eq]
			pair_flag_extract = Glob_pair_flag_cluster_temp[where_tray_eq]

			where_Xtop = np.where(Si_id_extract == 0)
			where_Xtop = where_Xtop[0]

			if len(where_Xtop) != 0:
				Si_id_intray = Si_id_extract[where_Xtop]
				tray_id_intray = tray_id_extract[where_Xtop]
				plane_id_intray = plane_id_extract[where_Xtop]
				pos_intray = pos_extract[where_Xtop]
				zpos_intray = zpos_extract[where_Xtop]
				energy_dep_intray = energy_dep_extract[where_Xtop]
				strip_number_intray = strip_number_extract[where_Xtop]
				pair_flag_intray = pair_flag_extract[where_Xtop]

				Si_id_temp = np.append(Si_id_temp, Si_id_intray)
				tray_id_temp = np.append(tray_id_temp, tray_id_intray)
				plane_id_temp = np.append(plane_id_temp, plane_id_intray)
				pos_temp = np.append(pos_temp, pos_intray)
				zpos_temp = np.append(zpos_temp, zpos_intray)
				energy_dep_temp = np.append(energy_dep_temp, energy_dep_intray)
				strip_number_temp = np.append(strip_number_temp, strip_number_intray)
				pair_flag_temp = np.append(pair_flag_temp, pair_flag_intray)				
					
			where_Ytop_arr = np.where(Si_id_extract == 1)
			where_Ytop = where_Ytop_arr[0]
					
					
			if len(where_Ytop) != 0:
				Si_id_intray = Si_id_extract[where_Ytop]
				tray_id_intray = tray_id_extract[where_Ytop]
				plane_id_intray = plane_id_extract[where_Ytop]
				pos_intray = pos_extract[where_Ytop]
				zpos_intray = zpos_extract[where_Ytop]
				energy_dep_intray = energy_dep_extract[where_Ytop]
				strip_number_intray = strip_number_extract[where_Ytop]
				pair_flag_intray = pair_flag_extract[where_Ytop]

				Si_id_temp = np.append(Si_id_temp, Si_id_intray)
				tray_id_temp = np.append(tray_id_temp, tray_id_intray)
				plane_id_temp = np.append(plane_id_temp, plane_id_intray)
				pos_temp = np.append(pos_temp, pos_intray)
				zpos_temp = np.append(zpos_temp, zpos_intray)
				energy_dep_temp = np.append(energy_dep_temp, energy_dep_intray)
				strip_number_temp = np.append(strip_number_temp, strip_number_intray)
				pair_flag_temp = np.append(pair_flag_temp, pair_flag_intray)				
					

			N_tray_eq = len(where_tray_eq)
			if where_tray_eq[N_tray_eq-1] < len(Glob_tray_id_cluster_temp)-1:
				intray = where_tray_eq[N_tray_eq-1]+1
			else:
				break

		event_id_temp = np.zeros(len(Si_id_temp), dtype=np.int64)
		for k in range(len(Si_id_temp)):
			event_id_temp[k] = event_array[j]

		Glob_event_id_cluster = np.append(Glob_event_id_cluster, event_id_temp)
		Glob_Si_id_cluster = np.append(Glob_Si_id_cluster, Si_id_temp)
		Glob_tray_id_cluster = np.append(Glob_tray_id_cluster, tray_id_temp)
		Glob_plane_id_cluster = np.append(Glob_plane_id_cluster, plane_id_temp)
		Glob_pos_cluster = np.append(Glob_pos_cluster, pos_temp)
		Glob_zpos_cluster = np.append(Glob_zpos_cluster, zpos_temp)
		Glob_energy_dep_cluster = np.append(Glob_energy_dep_cluster, energy_dep_temp)
		Glob_Strip_number_cluster = np.append(Glob_Strip_number_cluster, strip_number_temp)
		Glob_pair_flag_cluster = np.append(Glob_pair_flag_cluster, pair_flag_temp)

		j = j + 1

	Glob_event_id_cluster = np.array((Glob_event_id_cluster), dtype=np.int64)
	Glob_Si_id_cluster = np.array((Glob_Si_id_cluster), dtype=np.int64)
	Glob_tray_id_cluster = np.array((Glob_tray_id_cluster), dtype=np.int64)
	Glob_plane_id_cluster = np.array((Glob_plane_id_cluster), dtype=np.int64)
	Glob_pos_cluster = np.array(Glob_pos_cluster)
	Glob_zpos_cluster = np.array(Glob_zpos_cluster)
	Glob_energy_dep_cluster = np.array(Glob_energy_dep_cluster)
	Glob_Strip_number_cluster = np.array((Glob_Strip_number_cluster), dtype=np.int64)
	Glob_pair_flag_cluster = np.array((Glob_pair_flag_cluster), dtype=np.int64)

	return(Glob_event_id_cluster, Glob_Si_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_Strip_number_cluster, Glob_pair_flag_cluster)
