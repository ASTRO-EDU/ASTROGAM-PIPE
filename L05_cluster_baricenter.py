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


############# L0.5 cluster baricenter #################################################################################################################################################################

			
def L05_cluster_x(N_trig, Glob_plane_id_x_top, Glob_vol_id_x_top, Glob_moth_id_x_top, Glob_Strip_id_x_top, Glob_Si_id_x_top, Glob_tray_id_x_top, Glob_xpos_x_top, Glob_zpos_x_top, Glob_energy_dep_x_top, Glob_pair_flag_x_top): 		
			
	Glob_event_id_x_top_cluster = []
	Glob_Si_id_x_top_cluster = []
	Glob_tray_id_x_top_cluster = []
	Glob_plane_id_x_top_cluster = []
	Glob_zpos_x_top_cluster = []
	Glob_energy_dep_x_top_cluster = []
	Glob_xpos_x_top_cluster = []
	Glob_Strip_number_x_top_cluster = []
	Glob_pair_flag_x_top_cluster = []
			

	for k in range(N_trig):

		N_start = 0
		j=0
				
		while 1:

			############ Funzione di ordinamento 3
				
			sort_ascending_plane_x = np.argsort(Glob_plane_id_x_top[:,k], kind='mergesort')				
				
			############ fine funzione di ordinamento	


			Glob_vol_id_x_top_tray = Glob_vol_id_x_top[sort_ascending_plane_x, k]
			Glob_moth_id_x_top_tray = Glob_moth_id_x_top[sort_ascending_plane_x, k]
			Glob_Strip_id_x_top_tray = Glob_Strip_id_x_top[sort_ascending_plane_x, k]
			Glob_Si_id_x_top_tray = Glob_Si_id_x_top[sort_ascending_plane_x, k]
			Glob_tray_id_x_top_tray = Glob_tray_id_x_top[sort_ascending_plane_x, k]
			Glob_plane_id_x_top_tray = Glob_plane_id_x_top[sort_ascending_plane_x, k]
			Glob_xpos_x_top_tray = Glob_xpos_x_top[sort_ascending_plane_x, k]
			Glob_zpos_x_top_tray = Glob_zpos_x_top[sort_ascending_plane_x, k]
			Glob_energy_dep_x_top_tray = Glob_energy_dep_x_top[sort_ascending_plane_x, k]
			Glob_pair_flag_x_top_tray = Glob_pair_flag_x_top[sort_ascending_plane_x, k]

			where_tray_eq_x_top = np.where(Glob_tray_id_x_top_tray == Glob_tray_id_x_top_tray[j])
			where_tray_eq_x_top = where_tray_eq_x_top[0] 

			Glob_vol_id_x_top_tray = Glob_vol_id_x_top_tray[where_tray_eq_x_top]
			Glob_moth_id_x_top_tray = Glob_moth_id_x_top_tray[where_tray_eq_x_top]
			Glob_Strip_id_x_top_tray = Glob_Strip_id_x_top_tray[where_tray_eq_x_top]
			Glob_Si_id_x_top_tray = Glob_Si_id_x_top_tray[where_tray_eq_x_top]
			Glob_tray_id_x_top_tray = Glob_tray_id_x_top_tray[where_tray_eq_x_top]
			Glob_plane_id_x_top_tray = Glob_plane_id_x_top_tray[where_tray_eq_x_top]
			Glob_xpos_x_top_tray = Glob_xpos_x_top_tray[where_tray_eq_x_top]
			Glob_zpos_x_top_tray = Glob_zpos_x_top_tray[where_tray_eq_x_top]
			Glob_energy_dep_x_top_tray = Glob_energy_dep_x_top_tray[where_tray_eq_x_top]
			Glob_pair_flag_x_top_tray = Glob_pair_flag_x_top_tray[where_tray_eq_x_top]


			where_layer_x_top = np.where((Glob_Si_id_x_top_tray == 0) & (Glob_energy_dep_x_top_tray > 0.))
			where_layer_x_top = where_layer_x_top[0]


			if len(where_layer_x_top) != 0:
				Glob_vol_id_x_top_tray = Glob_vol_id_x_top_tray[where_layer_x_top]
				Glob_moth_id_x_top_tray = Glob_moth_id_x_top_tray[where_layer_x_top]
				Glob_Strip_id_x_top_tray = Glob_Strip_id_x_top_tray[where_layer_x_top]
				Glob_Si_id_x_top_tray = Glob_Si_id_x_top_tray[where_layer_x_top]
				Glob_tray_id_x_top_tray = Glob_tray_id_x_top_tray[where_layer_x_top]
				Glob_plane_id_x_top_tray = Glob_plane_id_x_top_tray[where_layer_x_top]
				Glob_xpos_x_top_tray = Glob_xpos_x_top_tray[where_layer_x_top]
				Glob_zpos_x_top_tray = Glob_zpos_x_top_tray[where_layer_x_top]
				Glob_energy_dep_x_top_tray = Glob_energy_dep_x_top_tray[where_layer_x_top]
				Glob_pair_flag_x_top_tray = Glob_pair_flag_x_top_tray[where_layer_x_top]


				############ Funzione di ordinamento 4

				sort_strip_ascending = np.argsort(Glob_Strip_id_x_top_tray, kind='mergesort')
									
				############ fine funzione di ordinamento	

	
				Glob_vol_id_x_top_tray = Glob_vol_id_x_top_tray[sort_strip_ascending]
				Glob_moth_id_x_top_tray = Glob_moth_id_x_top_tray[sort_strip_ascending]
				Glob_Strip_id_x_top_tray = Glob_Strip_id_x_top_tray[sort_strip_ascending]
				Glob_Si_id_x_top_tray = Glob_Si_id_x_top_tray[sort_strip_ascending]
				Glob_tray_id_x_top_tray = Glob_tray_id_x_top_tray[sort_strip_ascending]
				Glob_plane_id_x_top_tray = Glob_plane_id_x_top_tray[sort_strip_ascending]
				Glob_xpos_x_top_tray = Glob_xpos_x_top_tray[sort_strip_ascending]
				Glob_zpos_x_top_tray = Glob_zpos_x_top_tray[sort_strip_ascending]
				Glob_energy_dep_x_top_tray = Glob_energy_dep_x_top_tray[sort_strip_ascending]
				Glob_pair_flag_x_top_tray = Glob_pair_flag_x_top_tray[sort_strip_ascending]
			

	
				e_cluster_temp = [Glob_energy_dep_x_top_tray[0]]
				wx_cluster_temp = [Glob_xpos_x_top_tray[0]*Glob_energy_dep_x_top_tray[0]]
				nstrip_temp = 1


				if len(Glob_Strip_id_x_top_tray) == 1:

					Glob_event_id_x_top_cluster_old = k
					Glob_Si_id_x_top_cluster_old = Glob_Si_id_x_top_tray
					Glob_tray_id_x_top_cluster_old = Glob_tray_id_x_top_tray
					Glob_plane_id_x_top_cluster_old = Glob_plane_id_x_top_tray
					Glob_zpos_x_top_cluster_old = Glob_zpos_x_top_tray
					Glob_energy_dep_x_top_cluster_old = np.sum(e_cluster_temp)
					Glob_xpos_x_top_cluster_old = (np.sum(wx_cluster_temp)/np.sum(e_cluster_temp))
					Glob_Strip_number_x_top_cluster_old = nstrip_temp
					Glob_pair_flag_x_top_cluster_old = Glob_pair_flag_x_top_tray


					Glob_event_id_x_top_cluster = np.append(Glob_event_id_x_top_cluster, Glob_event_id_x_top_cluster_old)
					Glob_Si_id_x_top_cluster = np.append(Glob_Si_id_x_top_cluster, Glob_Si_id_x_top_cluster_old)
					Glob_tray_id_x_top_cluster = np.append(Glob_tray_id_x_top_cluster, Glob_tray_id_x_top_cluster_old)
					Glob_plane_id_x_top_cluster = np.append(Glob_plane_id_x_top_cluster, Glob_plane_id_x_top_cluster_old)
					Glob_zpos_x_top_cluster = np.append(Glob_zpos_x_top_cluster, Glob_zpos_x_top_cluster_old)
					Glob_energy_dep_x_top_cluster = np.append(Glob_energy_dep_x_top_cluster, Glob_energy_dep_x_top_cluster_old)
					Glob_xpos_x_top_cluster = np.append(Glob_xpos_x_top_cluster, Glob_xpos_x_top_cluster_old)
					Glob_Strip_number_x_top_cluster = np.append(Glob_Strip_number_x_top_cluster, Glob_Strip_number_x_top_cluster_old)
					Glob_pair_flag_x_top_cluster = np.append(Glob_pair_flag_x_top_cluster, Glob_pair_flag_x_top_cluster_old)

				else:

					for jc in range(len(Glob_Strip_id_x_top_tray)-1):

						if Glob_Strip_id_x_top_tray[jc+1] == (Glob_Strip_id_x_top_tray[jc]+1):
							e_cluster_temp_old = Glob_energy_dep_x_top_tray[jc+1]
							wx_cluster_temp_old = Glob_xpos_x_top_tray[jc + 1]*Glob_energy_dep_x_top_tray[jc+1]
							nstrip_temp = nstrip_temp+1						

							e_cluster_temp = np.append(e_cluster_temp, e_cluster_temp_old)
							wx_cluster_temp = np.append(wx_cluster_temp, wx_cluster_temp_old)


							if jc == len(Glob_Strip_id_x_top_tray)-2:
																			
								Glob_event_id_x_top_cluster_old = k
								Glob_Si_id_x_top_cluster_old = Glob_Si_id_x_top_tray[jc]
								Glob_tray_id_x_top_cluster_old = Glob_tray_id_x_top_tray[jc]
								Glob_plane_id_x_top_cluster_old = Glob_plane_id_x_top_tray[jc]
								Glob_zpos_x_top_cluster_old = Glob_zpos_x_top_tray[jc]
								Glob_energy_dep_x_top_cluster_old = np.sum(e_cluster_temp)
								Glob_xpos_x_top_cluster_old = (np.sum(wx_cluster_temp)/np.sum(e_cluster_temp))
								Glob_Strip_number_x_top_cluster_old = nstrip_temp
								Glob_pair_flag_x_top_cluster_old = Glob_pair_flag_x_top_tray[jc]


								Glob_event_id_x_top_cluster = np.append(Glob_event_id_x_top_cluster, Glob_event_id_x_top_cluster_old)
								Glob_Si_id_x_top_cluster = np.append(Glob_Si_id_x_top_cluster, Glob_Si_id_x_top_cluster_old)
								Glob_tray_id_x_top_cluster = np.append(Glob_tray_id_x_top_cluster, Glob_tray_id_x_top_cluster_old)
								Glob_plane_id_x_top_cluster = np.append(Glob_plane_id_x_top_cluster, Glob_plane_id_x_top_cluster_old)
								Glob_zpos_x_top_cluster = np.append(Glob_zpos_x_top_cluster, Glob_zpos_x_top_cluster_old)
								Glob_energy_dep_x_top_cluster = np.append(Glob_energy_dep_x_top_cluster, Glob_energy_dep_x_top_cluster_old)
								Glob_xpos_x_top_cluster = np.append(Glob_xpos_x_top_cluster, Glob_xpos_x_top_cluster_old)
								Glob_Strip_number_x_top_cluster = np.append(Glob_Strip_number_x_top_cluster, Glob_Strip_number_x_top_cluster_old)
								Glob_pair_flag_x_top_cluster = np.append(Glob_pair_flag_x_top_cluster, Glob_pair_flag_x_top_cluster_old)

						else:

							Glob_event_id_x_top_cluster_old = k
							Glob_Si_id_x_top_cluster_old = Glob_Si_id_x_top_tray[jc]
							Glob_tray_id_x_top_cluster_old = Glob_tray_id_x_top_tray[jc]
							Glob_plane_id_x_top_cluster_old = Glob_plane_id_x_top_tray[jc]
							Glob_zpos_x_top_cluster_old = Glob_zpos_x_top_tray[jc]
							Glob_energy_dep_x_top_cluster_old = np.sum(e_cluster_temp)
							Glob_xpos_x_top_cluster_old = (np.sum(wx_cluster_temp)/np.sum(e_cluster_temp))
							Glob_Strip_number_x_top_cluster_old = nstrip_temp
							Glob_pair_flag_x_top_cluster_old = Glob_pair_flag_x_top_tray[jc]
						
							Glob_event_id_x_top_cluster = np.append(Glob_event_id_x_top_cluster, Glob_event_id_x_top_cluster_old)
							Glob_Si_id_x_top_cluster = np.append(Glob_Si_id_x_top_cluster, Glob_Si_id_x_top_cluster_old)
							Glob_tray_id_x_top_cluster = np.append(Glob_tray_id_x_top_cluster, Glob_tray_id_x_top_cluster_old)
							Glob_plane_id_x_top_cluster = np.append(Glob_plane_id_x_top_cluster, Glob_plane_id_x_top_cluster_old)
							Glob_zpos_x_top_cluster = np.append(Glob_zpos_x_top_cluster, Glob_zpos_x_top_cluster_old)
							Glob_energy_dep_x_top_cluster = np.append(Glob_energy_dep_x_top_cluster, Glob_energy_dep_x_top_cluster_old)
							Glob_xpos_x_top_cluster = np.append(Glob_xpos_x_top_cluster, Glob_xpos_x_top_cluster_old)
							Glob_Strip_number_x_top_cluster = np.append(Glob_Strip_number_x_top_cluster, Glob_Strip_number_x_top_cluster_old)
							Glob_pair_flag_x_top_cluster = np.append(Glob_pair_flag_x_top_cluster, Glob_pair_flag_x_top_cluster_old)


							e_cluster_temp = [Glob_energy_dep_x_top_tray[jc+1]]
							wx_cluster_temp = [Glob_xpos_x_top_tray[jc+1]*Glob_energy_dep_x_top_tray[jc+1]]
							nstrip_temp = 1

							if jc == len(Glob_Strip_id_x_top_tray)-2:
								Glob_event_id_x_top_cluster_old = k
								Glob_Si_id_x_top_cluster_old = Glob_Si_id_x_top_tray[jc+1]
								Glob_tray_id_x_top_cluster_old = Glob_tray_id_x_top_tray[jc+1]
								Glob_plane_id_x_top_cluster_old = Glob_plane_id_x_top_tray[jc+1]
								Glob_zpos_x_top_cluster_old = Glob_zpos_x_top_tray[jc+1]
								Glob_energy_dep_x_top_cluster_old = Glob_energy_dep_x_top_tray[jc+1]
								Glob_xpos_x_top_cluster_old = Glob_xpos_x_top_tray[jc +1]
								Glob_Strip_number_x_top_cluster_old = nstrip_temp
								Glob_pair_flag_x_top_cluster_old = Glob_pair_flag_x_top_tray[jc+1]

								Glob_event_id_x_top_cluster = np.append(Glob_event_id_x_top_cluster, Glob_event_id_x_top_cluster_old)
								Glob_Si_id_x_top_cluster = np.append(Glob_Si_id_x_top_cluster, Glob_Si_id_x_top_cluster_old)
								Glob_tray_id_x_top_cluster = np.append(Glob_tray_id_x_top_cluster, Glob_tray_id_x_top_cluster_old)
								Glob_plane_id_x_top_cluster = np.append(Glob_plane_id_x_top_cluster, Glob_plane_id_x_top_cluster_old)
								Glob_zpos_x_top_cluster = np.append(Glob_zpos_x_top_cluster, Glob_zpos_x_top_cluster_old)
								Glob_energy_dep_x_top_cluster = np.append(Glob_energy_dep_x_top_cluster, Glob_energy_dep_x_top_cluster_old)
								Glob_xpos_x_top_cluster = np.append(Glob_xpos_x_top_cluster, Glob_xpos_x_top_cluster_old)
								Glob_Strip_number_x_top_cluster = np.append(Glob_Strip_number_x_top_cluster, Glob_Strip_number_x_top_cluster_old)
								Glob_pair_flag_x_top_cluster = np.append(Glob_pair_flag_x_top_cluster, Glob_pair_flag_x_top_cluster_old)


			N_tray_eq_x = len(where_tray_eq_x_top)
			if where_tray_eq_x_top[N_tray_eq_x-1] < np.size(Glob_tray_id_x_top[:,k])-1: # verificare con il len 
				j = where_tray_eq_x_top[N_tray_eq_x-1]+1
			else:
				break

	Glob_event_id_x_top_cluster = np.array((Glob_event_id_x_top_cluster), dtype=np.int64)
	Glob_Si_id_x_top_cluster = np.array((Glob_Si_id_x_top_cluster), dtype=np.int64)			
	Glob_tray_id_x_top_cluster = np.array((Glob_tray_id_x_top_cluster), dtype=np.int64)
	Glob_plane_id_x_top_cluster = np.array((Glob_plane_id_x_top_cluster), dtype=np.int64)
	Glob_zpos_x_top_cluster = np.array(Glob_zpos_x_top_cluster)
	Glob_energy_dep_x_top_cluster = np.array(Glob_energy_dep_x_top_cluster)
	Glob_xpos_x_top_cluster = np.array(Glob_xpos_x_top_cluster)
	Glob_Strip_number_x_top_cluster = np.array((Glob_Strip_number_x_top_cluster), dtype=np.int64)
	Glob_pair_flag_x_top_cluster = np.array(Glob_pair_flag_x_top_cluster)

	return(Glob_event_id_x_top_cluster, Glob_Si_id_x_top_cluster, Glob_tray_id_x_top_cluster, Glob_plane_id_x_top_cluster, Glob_zpos_x_top_cluster, Glob_energy_dep_x_top_cluster, Glob_xpos_x_top_cluster, Glob_Strip_number_x_top_cluster, Glob_pair_flag_x_top_cluster)
	
	
	

def L05_cluster_y(N_trig, Glob_plane_id_y_top, Glob_vol_id_y_top, Glob_moth_id_y_top, Glob_Strip_id_y_top, Glob_Si_id_y_top, Glob_tray_id_y_top, Glob_ypos_y_top, Glob_zpos_y_top, Glob_energy_dep_y_top, Glob_pair_flag_y_top): 		
	


	Glob_event_id_y_top_cluster = []
	Glob_Si_id_y_top_cluster = []
	Glob_tray_id_y_top_cluster = []
	Glob_plane_id_y_top_cluster = []
	Glob_zpos_y_top_cluster = []
	Glob_energy_dep_y_top_cluster = []
	Glob_ypos_y_top_cluster = []
	Glob_Strip_number_y_top_cluster = []
	Glob_pair_flag_y_top_cluster = []

			
	for k in range(N_trig):

		N_start = 0
		j=0

		while 1:


			############ Funzione di ordinamento 5
				
			sort_ascending_plane_y = np.argsort(Glob_plane_id_y_top[:,k], kind='mergesort')									
				
			############ fine funzione di ordinamento	


			Glob_vol_id_y_top_tray = Glob_vol_id_y_top[sort_ascending_plane_y, k]
			Glob_moth_id_y_top_tray = Glob_moth_id_y_top[sort_ascending_plane_y, k]
			Glob_Strip_id_y_top_tray = Glob_Strip_id_y_top[sort_ascending_plane_y, k]
			Glob_Si_id_y_top_tray = Glob_Si_id_y_top[sort_ascending_plane_y, k]
			Glob_tray_id_y_top_tray = Glob_tray_id_y_top[sort_ascending_plane_y, k]
			Glob_plane_id_y_top_tray = Glob_plane_id_y_top[sort_ascending_plane_y, k]
			Glob_ypos_y_top_tray = Glob_ypos_y_top[sort_ascending_plane_y, k]
			Glob_zpos_y_top_tray = Glob_zpos_y_top[sort_ascending_plane_y, k]
			Glob_energy_dep_y_top_tray = Glob_energy_dep_y_top[sort_ascending_plane_y, k]
			Glob_pair_flag_y_top_tray = Glob_pair_flag_y_top[sort_ascending_plane_y, k]

			where_tray_eq_y_top_arr = np.where(Glob_tray_id_y_top_tray == Glob_tray_id_y_top_tray[j])
			where_tray_eq_y_top = where_tray_eq_y_top_arr[0] 

			Glob_vol_id_y_top_tray = Glob_vol_id_y_top_tray[where_tray_eq_y_top]
			Glob_moth_id_y_top_tray = Glob_moth_id_y_top_tray[where_tray_eq_y_top]
			Glob_Strip_id_y_top_tray = Glob_Strip_id_y_top_tray[where_tray_eq_y_top]
			Glob_Si_id_y_top_tray = Glob_Si_id_y_top_tray[where_tray_eq_y_top]
			Glob_tray_id_y_top_tray = Glob_tray_id_y_top_tray[where_tray_eq_y_top]
			Glob_plane_id_y_top_tray = Glob_plane_id_y_top_tray[where_tray_eq_y_top]
			Glob_ypos_y_top_tray = Glob_ypos_y_top_tray[where_tray_eq_y_top]
			Glob_zpos_y_top_tray = Glob_zpos_y_top_tray[where_tray_eq_y_top]
			Glob_energy_dep_y_top_tray = Glob_energy_dep_y_top_tray[where_tray_eq_y_top]
			Glob_pair_flag_y_top_tray = Glob_pair_flag_y_top_tray[where_tray_eq_y_top]
					
					
			where_layer_y_top = np.where((Glob_Si_id_y_top_tray == 1) & (Glob_energy_dep_y_top_tray > 0.))
			where_layer_y_top = where_layer_y_top[0]
					
			if len(where_layer_y_top) != 0:
				Glob_vol_id_y_top_tray = Glob_vol_id_y_top_tray[where_layer_y_top]
				Glob_moth_id_y_top_tray = Glob_moth_id_y_top_tray[where_layer_y_top]
				Glob_Strip_id_y_top_tray = Glob_Strip_id_y_top_tray[where_layer_y_top]
				Glob_Si_id_y_top_tray = Glob_Si_id_y_top_tray[where_layer_y_top]
				Glob_tray_id_y_top_tray = Glob_tray_id_y_top_tray[where_layer_y_top]
				Glob_plane_id_y_top_tray = Glob_plane_id_y_top_tray[where_layer_y_top]
				Glob_ypos_y_top_tray = Glob_ypos_y_top_tray[where_layer_y_top]
				Glob_zpos_y_top_tray = Glob_zpos_y_top_tray[where_layer_y_top]
				Glob_energy_dep_y_top_tray = Glob_energy_dep_y_top_tray[where_layer_y_top]
				Glob_pair_flag_y_top_tray = Glob_pair_flag_y_top_tray[where_layer_y_top]
						
						
				############ Funzione di ordinamento 4
					
				sort_strip_ascending = np.argsort(Glob_Strip_id_y_top_tray, kind='mergesort')
				
				############ fine funzione di ordinamento	


				Glob_vol_id_y_top_tray = Glob_vol_id_y_top_tray[sort_strip_ascending]
				Glob_moth_id_y_top_tray = Glob_moth_id_y_top_tray[sort_strip_ascending]
				Glob_Strip_id_y_top_tray = Glob_Strip_id_y_top_tray[sort_strip_ascending]
				Glob_Si_id_y_top_tray = Glob_Si_id_y_top_tray[sort_strip_ascending]
				Glob_tray_id_y_top_tray = Glob_tray_id_y_top_tray[sort_strip_ascending]
				Glob_plane_id_y_top_tray = Glob_plane_id_y_top_tray[sort_strip_ascending]
				Glob_ypos_y_top_tray = Glob_ypos_y_top_tray[sort_strip_ascending]
				Glob_zpos_y_top_tray = Glob_zpos_y_top_tray[sort_strip_ascending]
				Glob_energy_dep_y_top_tray = Glob_energy_dep_y_top_tray[sort_strip_ascending]
				Glob_pair_flag_y_top_tray = Glob_pair_flag_y_top_tray[sort_strip_ascending]

	
				e_cluster_temp = [Glob_energy_dep_y_top_tray[0]]
				wx_cluster_temp = [Glob_ypos_y_top_tray[0]*Glob_energy_dep_y_top_tray[0]]
				nstrip_temp = 1


				if len(Glob_Strip_id_y_top_tray) == 1:

					Glob_event_id_y_top_cluster_old = k
					Glob_Si_id_y_top_cluster_old = Glob_Si_id_y_top_tray
					Glob_tray_id_y_top_cluster_old = Glob_tray_id_y_top_tray
					Glob_plane_id_y_top_cluster_old = Glob_plane_id_y_top_tray
					Glob_zpos_y_top_cluster_old = Glob_zpos_y_top_tray
					Glob_energy_dep_y_top_cluster_old = np.sum(e_cluster_temp)
					Glob_ypos_y_top_cluster_old = (np.sum(wx_cluster_temp)/np.sum(e_cluster_temp))
					Glob_Strip_number_y_top_cluster_old = nstrip_temp
					Glob_pair_flag_y_top_cluster_old = Glob_pair_flag_y_top_tray

					Glob_event_id_y_top_cluster = np.append(Glob_event_id_y_top_cluster, Glob_event_id_y_top_cluster_old)
					Glob_Si_id_y_top_cluster = np.append(Glob_Si_id_y_top_cluster, Glob_Si_id_y_top_cluster_old)
					Glob_tray_id_y_top_cluster = np.append(Glob_tray_id_y_top_cluster, Glob_tray_id_y_top_cluster_old)
					Glob_plane_id_y_top_cluster = np.append(Glob_plane_id_y_top_cluster, Glob_plane_id_y_top_cluster_old)
					Glob_zpos_y_top_cluster = np.append(Glob_zpos_y_top_cluster, Glob_zpos_y_top_cluster_old)
					Glob_energy_dep_y_top_cluster = np.append(Glob_energy_dep_y_top_cluster, Glob_energy_dep_y_top_cluster_old)
					Glob_ypos_y_top_cluster = np.append(Glob_ypos_y_top_cluster, Glob_ypos_y_top_cluster_old)
					Glob_Strip_number_y_top_cluster = np.append(Glob_Strip_number_y_top_cluster, Glob_Strip_number_y_top_cluster_old)
					Glob_pair_flag_y_top_cluster = np.append(Glob_pair_flag_y_top_cluster, Glob_pair_flag_y_top_cluster_old)

				else:
				
					for jc in range(len(Glob_Strip_id_y_top_tray)-1):

						if Glob_Strip_id_y_top_tray[jc+1] == (Glob_Strip_id_y_top_tray[jc]+1):
							e_cluster_temp_old = Glob_energy_dep_y_top_tray[jc+1]
							wx_cluster_temp_old = Glob_ypos_y_top_tray[jc + 1]*Glob_energy_dep_y_top_tray[jc+1]
							nstrip_temp = nstrip_temp+1						

							e_cluster_temp = np.append(e_cluster_temp, e_cluster_temp_old)
							wx_cluster_temp = np.append(wx_cluster_temp, wx_cluster_temp_old)

								
							if jc == len(Glob_Strip_id_y_top_tray)-2:
																		
								Glob_event_id_y_top_cluster_old = k
								Glob_Si_id_y_top_cluster_old = Glob_Si_id_y_top_tray[jc]
								Glob_tray_id_y_top_cluster_old = Glob_tray_id_y_top_tray[jc]
								Glob_plane_id_y_top_cluster_old = Glob_plane_id_y_top_tray[jc]
								Glob_zpos_y_top_cluster_old = Glob_zpos_y_top_tray[jc]
								Glob_energy_dep_y_top_cluster_old = np.sum(e_cluster_temp)
								Glob_ypos_y_top_cluster_old = (np.sum(wx_cluster_temp)/np.sum(e_cluster_temp))
								Glob_Strip_number_y_top_cluster_old = nstrip_temp
								Glob_pair_flag_y_top_cluster_old = Glob_pair_flag_y_top_tray[jc]

								Glob_event_id_y_top_cluster = np.append(Glob_event_id_y_top_cluster, Glob_event_id_y_top_cluster_old)
								Glob_Si_id_y_top_cluster = np.append(Glob_Si_id_y_top_cluster, Glob_Si_id_y_top_cluster_old)
								Glob_tray_id_y_top_cluster = np.append(Glob_tray_id_y_top_cluster, Glob_tray_id_y_top_cluster_old)
								Glob_plane_id_y_top_cluster = np.append(Glob_plane_id_y_top_cluster, Glob_plane_id_y_top_cluster_old)
								Glob_zpos_y_top_cluster = np.append(Glob_zpos_y_top_cluster, Glob_zpos_y_top_cluster_old)
								Glob_energy_dep_y_top_cluster = np.append(Glob_energy_dep_y_top_cluster, Glob_energy_dep_y_top_cluster_old)
								Glob_ypos_y_top_cluster = np.append(Glob_ypos_y_top_cluster, Glob_ypos_y_top_cluster_old)
								Glob_Strip_number_y_top_cluster = np.append(Glob_Strip_number_y_top_cluster, Glob_Strip_number_y_top_cluster_old)
								Glob_pair_flag_y_top_cluster = np.append(Glob_pair_flag_y_top_cluster, Glob_pair_flag_y_top_cluster_old)

						else:

							Glob_event_id_y_top_cluster_old = k
							Glob_Si_id_y_top_cluster_old = Glob_Si_id_y_top_tray[jc]
							Glob_tray_id_y_top_cluster_old = Glob_tray_id_y_top_tray[jc]
							Glob_plane_id_y_top_cluster_old = Glob_plane_id_y_top_tray[jc]
							Glob_zpos_y_top_cluster_old = Glob_zpos_y_top_tray[jc]
							Glob_energy_dep_y_top_cluster_old = np.sum(e_cluster_temp)
							Glob_ypos_y_top_cluster_old = (np.sum(wx_cluster_temp)/np.sum(e_cluster_temp))
							Glob_Strip_number_y_top_cluster_old = nstrip_temp
							Glob_pair_flag_y_top_cluster_old = Glob_pair_flag_y_top_tray[jc]
									
							Glob_event_id_y_top_cluster = np.append(Glob_event_id_y_top_cluster, Glob_event_id_y_top_cluster_old)
							Glob_Si_id_y_top_cluster = np.append(Glob_Si_id_y_top_cluster, Glob_Si_id_y_top_cluster_old)
							Glob_tray_id_y_top_cluster = np.append(Glob_tray_id_y_top_cluster, Glob_tray_id_y_top_cluster_old)
							Glob_plane_id_y_top_cluster = np.append(Glob_plane_id_y_top_cluster, Glob_plane_id_y_top_cluster_old)
							Glob_zpos_y_top_cluster = np.append(Glob_zpos_y_top_cluster, Glob_zpos_y_top_cluster_old)
							Glob_energy_dep_y_top_cluster = np.append(Glob_energy_dep_y_top_cluster, Glob_energy_dep_y_top_cluster_old)
							Glob_ypos_y_top_cluster = np.append(Glob_ypos_y_top_cluster, Glob_ypos_y_top_cluster_old)
							Glob_Strip_number_y_top_cluster = np.append(Glob_Strip_number_y_top_cluster, Glob_Strip_number_y_top_cluster_old)
							Glob_pair_flag_y_top_cluster = np.append(Glob_pair_flag_y_top_cluster, Glob_pair_flag_y_top_cluster_old)


							e_cluster_temp = [Glob_energy_dep_y_top_tray[jc+1]]
							wx_cluster_temp = [Glob_ypos_y_top_tray[jc+1]*Glob_energy_dep_y_top_tray[jc+1]]
							nstrip_temp = 1

							if jc == len(Glob_Strip_id_y_top_tray)-2:
								Glob_event_id_y_top_cluster_old = k
								Glob_Si_id_y_top_cluster_old = Glob_Si_id_y_top_tray[jc+1]
								Glob_tray_id_y_top_cluster_old = Glob_tray_id_y_top_tray[jc+1]
								Glob_plane_id_y_top_cluster_old = Glob_plane_id_y_top_tray[jc+1]
								Glob_zpos_y_top_cluster_old = Glob_zpos_y_top_tray[jc+1]
								Glob_energy_dep_y_top_cluster_old = Glob_energy_dep_y_top_tray[jc+1]
								Glob_ypos_y_top_cluster_old = Glob_ypos_y_top_tray[jc +1]
								Glob_Strip_number_y_top_cluster_old = nstrip_temp
								Glob_pair_flag_y_top_cluster_old = Glob_pair_flag_y_top_tray[jc+1]

								Glob_event_id_y_top_cluster = np.append(Glob_event_id_y_top_cluster, Glob_event_id_y_top_cluster_old)
								Glob_Si_id_y_top_cluster = np.append(Glob_Si_id_y_top_cluster, Glob_Si_id_y_top_cluster_old)
								Glob_tray_id_y_top_cluster = np.append(Glob_tray_id_y_top_cluster, Glob_tray_id_y_top_cluster_old)
								Glob_plane_id_y_top_cluster = np.append(Glob_plane_id_y_top_cluster, Glob_plane_id_y_top_cluster_old)
								Glob_zpos_y_top_cluster = np.append(Glob_zpos_y_top_cluster, Glob_zpos_y_top_cluster_old)
								Glob_energy_dep_y_top_cluster = np.append(Glob_energy_dep_y_top_cluster, Glob_energy_dep_y_top_cluster_old)
								Glob_ypos_y_top_cluster = np.append(Glob_ypos_y_top_cluster, Glob_ypos_y_top_cluster_old)
								Glob_Strip_number_y_top_cluster = np.append(Glob_Strip_number_y_top_cluster, Glob_Strip_number_y_top_cluster_old)
								Glob_pair_flag_y_top_cluster = np.append(Glob_pair_flag_y_top_cluster, Glob_pair_flag_y_top_cluster_old)
										
										
			N_tray_eq_y = len(where_tray_eq_y_top)
			if where_tray_eq_y_top[N_tray_eq_y-1] < np.size(Glob_tray_id_y_top[:,k])-1: # verificare con il len 
				j = where_tray_eq_y_top[N_tray_eq_y-1]+1
			else:
				break

	Glob_event_id_y_top_cluster = np.array((Glob_event_id_y_top_cluster), dtype=np.int64)
	Glob_Si_id_y_top_cluster = np.array((Glob_Si_id_y_top_cluster), dtype=np.int64)			
	Glob_tray_id_y_top_cluster = np.array((Glob_tray_id_y_top_cluster), dtype=np.int64)
	Glob_plane_id_y_top_cluster = np.array((Glob_plane_id_y_top_cluster), dtype=np.int64)
	Glob_zpos_y_top_cluster = np.array(Glob_zpos_y_top_cluster)
	Glob_energy_dep_y_top_cluster = np.array(Glob_energy_dep_y_top_cluster)
	Glob_ypos_y_top_cluster = np.array(Glob_ypos_y_top_cluster)
	Glob_Strip_number_y_top_cluster = np.array((Glob_Strip_number_y_top_cluster), dtype=np.int64)
	Glob_pair_flag_y_top_cluster = np.array(Glob_pair_flag_y_top_cluster)


	return(Glob_event_id_y_top_cluster, Glob_Si_id_y_top_cluster, Glob_tray_id_y_top_cluster, Glob_plane_id_y_top_cluster, Glob_zpos_y_top_cluster, Glob_energy_dep_y_top_cluster, Glob_ypos_y_top_cluster, Glob_Strip_number_y_top_cluster, Glob_pair_flag_y_top_cluster)
	 	
###########################################################################################################################################################################################	

