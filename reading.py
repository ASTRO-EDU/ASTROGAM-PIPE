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

def reading_fits(filepath, ifile, cal_flag, ac_flag, tracker_top_vol_start, tracker_bottom_vol_start, tracker_top_bot_diff, cal_vol_start, cal_vol_end, ac_vol_start, ac_vol_end, part_type):

	t = fits.open(filepath+'/xyz.'+str(ifile)+'.fits.gz')
   
	tbdata = t[1].data

	evt_id = tbdata.field('EVT_ID')
	trk_id = tbdata.field('TRK_ID')
	parent_trk_id = tbdata.field('PARENT_TRK_ID')
	volume_id = tbdata.field('VOLUME_ID')
	volume_name = tbdata.field('VOLUME_NAME')
	mother_id = tbdata.field('MOTHER_ID')
	e_dep = tbdata.field('E_DEP')
	x_ent = tbdata.field('X_ENT')
	y_ent = tbdata.field('Y_ENT')
	z_ent = tbdata.field('Z_ENT')
	x_exit = tbdata.field('X_EXIT')
	y_exit = tbdata.field('Y_EXIT')
	z_exit = tbdata.field('Z_EXIT')
	e_kin_ent = tbdata.field('E_KIN_ENT')
	e_kin_exit = tbdata.field('E_KIN_EXIT')
	mdx_ent = tbdata.field('MDX_ENT')
	mdy_ent = tbdata.field('MDY_ENT')
	mdz_ent = tbdata.field('MDZ_ENT')
	mdx_exit = tbdata.field('MDX_EXIT')
	mdy_exit = tbdata.field('MDY_EXIT')
	mdz_exit = tbdata.field('MDZ_EXIT')
	gtime_ent = tbdata.field('GTIME_ENT')
	gtime_exit = tbdata.field('GTIME_EXIT')
	particle_id = tbdata.field('PARTICLE_ID')
	particle_name = tbdata.field('PARTICLE_NAME')
	process_id = tbdata.field('PROCESS_ID')
	process_name = tbdata.field('PROCESS_NAME')

	
	vol_id_tr = []	
	moth_id_tr = []
	event_id_tr = []
	en_dep_tr = []
	x_en_tr = []
	y_en_tr = []
	z_en_tr = []
	x_ex_tr = []
	y_ex_tr = []
	z_ex_tr = []
	child_id_tr = []
	proc_id_tr = []
	theta_ent_tr = []
	phi_ent_tr = []
	theta_exit_tr = []
	phi_exit_tr = []	
	x_pos = []
	y_pos = []
	z_pos = []

	event_id_cal = []
	moth_id_cal = []
	vol_id_cal = []
	energy_dep_cal = []
	x_en_cal = []
	y_en_cal = []
	z_en_cal = []
	x_ex_cal = []
	y_ex_cal = []
	z_ex_cal = []
	theta_ent_cal = []
	phi_ent_cal = []
	theta_exit_cal = []
	phi_exit_cal = []
	child_id_cal = []
	proc_id_cal = []			

	event_id_ac = []
	vol_id_ac = []
	en_dep_ac = []
	x_en_ac = []
	y_en_ac = []
	z_en_ac = []
	x_ex_ac = []
	y_ex_ac = []
	z_ex_ac = []
	theta_ent_ac = []
	phi_ent_ac = []
	theta_exit_ac = []
	phi_exit_ac = []
	child_id_ac = []
	proc_id_ac = []		
	moth_id_ac = []

	i = 0
	while i < len(tbdata):
		
		vol_id = volume_id[i]        
		moth_id = mother_id[i]
	        energy_dep = e_dep[i]         
	        cos_x_angle_ent = mdx_ent[i]
	        cos_y_angle_ent = mdy_ent[i]
	        cos_z_angle_ent = mdz_ent[i]
	        cos_x_angle_exit = mdx_exit[i]
	        cos_y_angle_exit = mdy_exit[i]
	        cos_z_angle_exit = mdz_exit[i]
	        gtime_en = gtime_ent[i]
	        gtime_ex = gtime_exit[i]
	        part_id = particle_id[i]
	        part_name = particle_name[i]
	        proc_name = process_name[i]
		
		#  Reading the tracker (events with E > 0)
       		        
		if vol_id >= tracker_bottom_vol_start or moth_id >= tracker_bottom_vol_start:

			if part_type == 'g':

				theta_ent = (180./math.pi)*math.acos(-(cos_x_angle_ent))
				phi_ent = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
	
				theta_exit = (180./math.pi)*math.acos(-(cos_z_angle_exit))
				phi_exit = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))
			
				x_en = (x_ent[i])/10.
				y_en = (y_ent[i])/10.
				z_en = (z_ent[i])/10.
				x_ex = (x_exit[i])/10.
				y_ex = (y_exit[i])/10.
				z_ex = (z_exit[i])/10.
			
				x_position = x_en + ((x_ex - x_en)/2.)
				y_position = y_en + ((y_ex - y_en)/2.)
				z_position = z_en + ((z_ex - z_en)/2.)

				vol_id_tr.append(volume_id[i])
				moth_id_tr.append(mother_id[i])
				event_id_tr.append(evt_id[i])
				en_dep_tr.append(100.)			
				x_en_tr.append((x_ent[i])/10.)
				y_en_tr.append((y_ent[i])/10.)
				z_en_tr.append((z_ent[i])/10.)
				x_ex_tr.append((x_exit[i])/10.)
				y_ex_tr.append((y_exit[i])/10.)
				z_ex_tr.append((z_exit[i])/10.)
				child_id_tr.append(parent_trk_id[i])
				proc_id_tr.append(process_id[i])
				theta_ent_tr.append(theta_ent)
				phi_ent_tr.append(phi_ent)
				theta_exit_tr.append(theta_exit)
				phi_exit_tr.append(phi_exit)
				x_pos.append(x_position)
				y_pos.append(y_position)
				z_pos.append(z_position)

			if energy_dep > 0.:

				theta_ent = (180./math.pi)*math.acos(-(cos_x_angle_ent))
				phi_ent = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
	
				theta_exit = (180./math.pi)*math.acos(-(cos_z_angle_exit))
				phi_exit = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

				x_en = (x_ent[i])/10.
				y_en = (y_ent[i])/10.
				z_en = (z_ent[i])/10.
				x_ex = (x_exit[i])/10.
				y_ex = (y_exit[i])/10.
				z_ex = (z_exit[i])/10.
			
				x_position = x_en + ((x_ex - x_en)/2.)
				y_position = y_en + ((y_ex - y_en)/2.)
				z_position = z_en + ((z_ex - z_en)/2.)

				vol_id_tr.append(volume_id[i])
				moth_id_tr.append(mother_id[i])
				event_id_tr.append(evt_id[i])
				en_dep_tr.append(e_dep[i])			
				x_en_tr.append((x_ent[i])/10.)
				y_en_tr.append((y_ent[i])/10.)
				z_en_tr.append((z_ent[i])/10.)
				x_ex_tr.append((x_exit[i])/10.)
				y_ex_tr.append((y_exit[i])/10.)
				z_ex_tr.append((z_exit[i])/10.)
				child_id_tr.append(parent_trk_id[i])
				proc_id_tr.append(process_id[i])
				theta_ent_tr.append(theta_ent)
				phi_ent_tr.append(phi_ent)
				theta_exit_tr.append(theta_exit)
				phi_exit_tr.append(phi_exit)
				x_pos.append(x_position)
				y_pos.append(y_position)
				z_pos.append(z_position)
				

		# Reading the Calorimeter
		if cal_flag == 1:

			if vol_id >= cal_vol_start and vol_id <= cal_vol_end:
				if part_type == 'g':

					theta_ent_calor = (180./math.pi)*math.acos(-(cos_x_angle_ent))
					phi_ent_calor = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
		
					theta_exit_calor = (180./math.pi)*math.acos(-(cos_z_angle_exit))
					phi_exit_calor = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

					event_id_cal.append(evt_id[i])
					moth_id_cal.append(mother_id[i])
					vol_id_cal.append(volume_id[i])
					energy_dep_cal.append(e_dep[i])
					x_en_cal.append((x_ent[i])/10.)
					y_en_cal.append((y_ent[i])/10.)
					z_en_cal.append((z_ent[i])/10.)
					x_ex_cal.append((x_exit[i])/10.)
					y_ex_cal.append((y_exit[i])/10.)
					z_ex_cal.append((z_exit[i])/10.)

					theta_ent_cal.append(theta_ent_calor)
					phi_ent_cal.append(phi_ent_calor)
		
					theta_exit_cal.append(theta_exit_calor)
					phi_exit_cal.append(phi_exit_calor)

					child_id_cal.append(parent_trk_id[i])
					proc_id_cal.append(process_id[i])
					

				if energy_dep > 0.:
					theta_ent_calor = (180./math.pi)*math.acos(-(cos_x_angle_ent))
					phi_ent_calor = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
		
					theta_exit_calor = (180./math.pi)*math.acos(-(cos_z_angle_exit))
					phi_exit_calor = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

					event_id_cal.append(evt_id[i])
					moth_id_cal.append(mother_id[i])
					vol_id_cal.append(volume_id[i])
					energy_dep_cal.append(e_dep[i])
					x_en_cal.append((x_ent[i])/10.)
					y_en_cal.append((y_ent[i])/10.)
					z_en_cal.append((z_ent[i])/10.)
					x_ex_cal.append((x_exit[i])/10.)
					y_ex_cal.append((y_exit[i])/10.)
					z_ex_cal.append((z_exit[i])/10.)

					theta_ent_cal.append(theta_ent_calor)
					phi_ent_cal.append(phi_ent_calor)
		
					theta_exit_cal.append(theta_exit_calor)
					phi_exit_cal.append(phi_exit_calor)

					child_id_cal.append(parent_trk_id[i])
					proc_id_cal.append(process_id[i])


		# Reading the AC
		if ac_flag == 1:

			if vol_id >= ac_vol_start and vol_id <= ac_vol_end:
	
				if part_type == 'g':

					theta_ent_anc = (180./math.pi)*math.acos(-(cos_x_angle_ent))
					phi_ent_anc = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
	
					theta_exit_anc = (180./math.pi)*math.acos(-(cos_z_angle_exit))
					phi_exit_anc = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

					if isStrip == 1:
						moth_id_ac.append(mother_id[i])
					else:
						moth_id_ac.append(0)

					event_id_ac.append(evt_id[i])
					vol_id_ac.append(vol_id[i])
					x_en_ac.append(x_ent[i]/10.)
					y_en_ac.append(y_ent[i]/10.)
					z_en_ac.append(z_ent[i]/10.)
					x_ex_ac.append(x_exit[i]/10.)
					y_ex_ac.append(y_exit[i]/10.)
					z_ex_ac.append(z_exit[i]/10.)
					en_dep_ac.append(e_dep[i])
					theta_ent_ac.append(theta_ent_anc)
					phi_ent_ac.append(phi_ent_anc)
		
					theta_exit_ac.append(theta_exit_anc)
					phi_exit_ac_.append(phi_exit_anc)

					child_id_ac.append(parent_trk_id[i])
					proc_id_ac.append(process_id[i])


				if energy_dep > 0.:
					
					theta_ent_anc = (180./math.pi)*math.acos(-(cos_x_angle_ent))
					phi_ent_anc = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
	
					theta_exit_anc = (180./math.pi)*math.acos(-(cos_z_angle_exit))
					phi_exit_anc = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

					if isStrip == 1:
						moth_id_ac.append(mother_id[i])
					else:
						moth_id_ac.append(0)

					event_id_ac.append(evt_id[i])
					vol_id_ac.append(vol_id[i])
					x_en_ac.append(x_ent[i]/10.)
					y_en_ac.append(y_ent[i]/10.)
					z_en_ac.append(z_ent[i]/10.)
					x_ex_ac.append(x_exit[i]/10.)
					y_ex_ac.append(y_exit[i]/10.)
					z_ex_ac.append(z_exit[i]/10.)
					en_dep_ac.append(e_dep[i])
					theta_ent_ac.append(theta_ent_anc)
					phi_ent_ac.append(phi_ent_anc)
		
					theta_exit_ac.append(theta_exit_anc)
					phi_exit_ac_.append(phi_exit_anc)

					child_id_ac.append(parent_trk_id[i])
					proc_id_ac.append(process_id[i])

	

		i = i + 1

	t.close()
	
###tracker arrays

	event_id_tr = np.array(event_id_tr)
	vol_id_tr = np.array(vol_id_tr)
	moth_id_tr = np.array(moth_id_tr)
	en_dep_tr = np.array(en_dep_tr)
	x_en_tr = np.array(x_en_tr)
	y_en_tr = np.array(y_en_tr)
	z_en_tr = np.array(z_en_tr)
	x_ex_tr = np.array(x_ex_tr)
	y_ex_tr = np.array(y_ex_tr)
	z_ex_tr = np.array(z_ex_tr)	
	child_id_tr = np.array(child_id_tr)
	proc_id_tr = np.array(proc_id_tr)
	x_pos = np.array(x_pos)
	y_pos = np.array(y_pos)
	z_pos = np.array(z_pos)

###### calorimeter arrays

	event_id_cal = np.array(event_id_cal)
	vol_id_cal = np.array(vol_id_cal)
	energy_dep_cal = np.array(energy_dep_cal)
	x_en_cal = np.array(x_en_cal)
	y_en_cal = np.array(y_en_cal)
	z_en_cal = np.array(z_en_cal)
	x_ex_cal = np.array(x_ex_cal)
	y_ex_cal = np.array(y_ex_cal)
	z_ex_cal = np.array(z_ex_cal)
	theta_ent_cal = np.array(theta_ent_cal)
	phi_ent_cal = np.array(phi_ent_cal)
	theta_exit_cal = np.array(theta_exit_cal)
	phi_exit_cal = np.array(phi_exit_cal)
	child_id_cal = np.array(child_id_cal)
	proc_id_cal = np.array(proc_id_cal)			
	moth_id_cal = np.array(moth_id_cal)

#### AC arrays

	event_id_ac = np.array(event_id_ac)
	vol_id_ac = np.array(vol_id_ac)
	en_dep_ac = np.array(en_dep_ac)
	moth_id_ac =np.array(moth_id_ac)


	return(vol_id_tr, moth_id_tr, event_id_tr, en_dep_tr, x_en_tr, y_en_tr, z_en_tr, x_ex_tr, y_ex_tr, z_ex_tr, child_id_tr, proc_id_tr, theta_ent_tr, phi_ent_tr, theta_exit_tr, phi_exit_tr, x_pos, y_pos, z_pos, vol_id_cal, moth_id_cal, event_id_cal, energy_dep_cal, x_en_cal, y_en_cal, z_en_cal, x_ex_cal, y_ex_cal, z_ex_cal, child_id_cal, proc_id_cal, theta_ent_cal, phi_ent_cal, theta_exit_cal, phi_exit_cal, vol_id_ac, event_id_ac, en_dep_ac, x_en_ac, y_en_ac, z_en_ac, x_ex_ac, y_ex_ac, z_ex_ac, child_id_ac, proc_id_ac, theta_ent_ac, phi_ent_ac, theta_exit_ac, phi_exit_ac)



	

