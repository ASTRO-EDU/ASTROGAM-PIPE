from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import pickle
from astropy.table import Table
import copy
import time



##############################
#     parametri iniziali     #
##############################

astrogam_version = 'V1.0'   # Enter eASTROGAM release (e.g. V1.0):
bogemms_tag = 211           # Enter BoGEMMS release (e.g. 211):
sim_type = 0                # Enter simulation type [0 = Mono, 1 = Range, 2 = Chen, 3: Vela, 4: Crab, 5: G400]:
py_list = 400               # Enter the Physics List [0 = QGSP_BERT_EMV, 100 = ARGO, 300 = FERMI, 400 = ASTROMEV]:
N_in = 100000               # Enter the number of emitted particles:
part_type = "ph"            # Enter the particle type [ph = photons, mu = muons, g = geantino, p = proton, el = electron]:
#n_fits = 1                  # Enter number of FITS files:
ene_range = 0               # Enter energy distribution [0 = MONO, 1 = POW, 2 = EXP, 3 = LIN]:
ene_min = 100               # Enter miminum energy [MeV]:
ene_max = 100               # Enter maximum energy [MeV]:
ang_type = "UNI"            # Enter the angular distribution [e.g. UNI, ISO]:
theta_type = 30             # Enter theta:
phi_type = 225              # Enter phi:
pol_type = 0                # Is the source polarized? [0 = false, 1 = true]:
pol_angle = 20              # Enter Polarization angle:
source_g = 0                # Enter source geometry [0 = Point, 1 = Plane]:
isStrip = 1                 # Strip/Pixels activated? [0 = false, 1 = true]
repli = 1                   # Strips/Pixels replicated? [0 = false, 1 = true]
cal_flag = 1                # Is Cal present? [0 = false, 1 = true]:
ac_flag = 0                 # Is AC present? [0 = false, 1 = true]:
passive_flag = 0            # Is Passive present? [0 = false, 1 = true]:
energy_thresh = 15          # Enter energy threshold [keV]:

if sim_type != 0 and sim_type != 1 and sim_type != 2 and sim_type != 3 and sim_type != 4 and sim_type != 5:
	exit('Error: sim_type could be 0 (Mono) - 1 (Range) - 2 (Chen) - 3 (Vela) - 4 (Crab) - 5 (G400)')

if py_list != 0 and py_list != 100 and py_list != 300 and py_list != 400:
	exit('Error: py_list could be 0 (QGSP_BERT_EMV) - 100 (ARGO) - 300 (FERMI) - 400 (ASTROMEV)')

if part_type != "ph" and part_type != "mu" and part_type != "g" and part_type != "p" and part_type != "el":
	exit('Error: part_type be ph (photons) - mu (muons) - g (geantino) - p (proton) - el (electron)')

if ene_range != 0 and ene_range != 1 and ene_range != 2 and ene_range != 3:
	exit('Error: ene_range could be 0 (MONO) - 1 (POW) - 2 (EXP) - 3 (LIN)')

if pol_type != 0 and pol_type != 1:
	exit('Error: pol_type could be 0 (false) - 1 (true)')

if source_g != 0 and source_g != 1:
	exit('Error: source_g could be 0 (Point) - 1 (Plane)')

if isStrip != 0 and isStrip != 1:
	exit('Error: isStrip could be 0 (false) - 1 (true)')

if repli != 0 and repli != 1:
	exit('Error: repli could be 0 (false) - 1 (true)')

if cal_flag != 0 and cal_flag != 1:
	exit('Error: cal_flag could be 0 (false) - 1 (true)')

if ac_flag != 0 and ac_flag != 1:
	exit('Error: ac_flag could be 0 (false) - 1 (true)')

if ac_flag != 0 and ac_flag != 1:
	exit('Error: ac_flag could be 0 (false) - 1 (true)')

if astrogam_version=='V1.0':
	astrogam_tag = '01'

sim_tag = 'eAST'+str(bogemms_tag)+str(astrogam_tag)+'0102'

if ene_range == 0:
	ene_dis = 'MONO'
	ene_type = ene_min
	if ene_type >= 1:
		ene_type = repr(ene_type)
	if ene_type < 1:
		ene_type = repr(ene_type)
		ene_type = ene_type[:5]
	if type(ene_type) is int:
		pass
	else:	
		nstring = len(ene_type)
		ene_type_notzero = ene_type
      		flag = 1 
		
		#for ichar_reverse in range(0, nstring):
			#ichar = (nstring-1) - ichar_reverse
		if ene_type[0] == '0' or  ene_type[0] == '.':
			if flag == 1:
				ene_type_notzero = ene_type_notzero[:5]
		else:
			flag = 0
		ene_type = ene_type_notzero
	


if ene_range == 1:
	ene_dis = 'POW'
	
	ene_min_string = repr(ene_min)	
	if type(ene_min) is int:
		pass
	else:
		nstring = len(ene_min_string)
		ene_min_string_notzero = ene_min_string
      		flag = 1 	
		if ene_min_string[0] == '0' or  ene_min_string[0] == '.':
			if flag == 1:
				ene_min_string_notzero = ene_min_string_notzero[:5]
		else:
			flag = 0
		ene_min_string = ene_min_string_notzero

	
	ene_max_string = repr(ene_max)	
	if type(ene_max) is int:
		pass
	else:
		nstring = len(ene_max_string)
		ene_max_string_notzero = ene_max_string
      		flag = 1 	
		if ene_max_string[0] == '0' or  ene_max_string[0] == '.':
			if flag == 1:
				ene_max_string_notzero = ene_max_string_notzero[:5]
		else:
			flag = 0
		ene_max_string = ene_max_string_notzero	
	
	ene_type = str(ene_min_string)+'-'+str(ene_max_string)
	



if ene_range == 2:
	ene_dis = 'EXP'
	
	ene_min_string = repr(ene_min)	
	if type(ene_min) is int:
		pass
	else:
		nstring = len(ene_min_string)
		ene_min_string_notzero = ene_min_string
      		flag = 1 	
		if ene_min_string[0] == '0' or  ene_min_string[0] == '.':
			if flag == 1:
				ene_min_string_notzero = ene_min_string_notzero[:5]
		else:
			flag = 0
		ene_min_string = ene_min_string_notzero

	
	ene_max_string = repr(ene_max)	
	if type(ene_max) is int:
		pass
	else:
		nstring = len(ene_max_string)
		ene_max_string_notzero = ene_max_string
      		flag = 1 	
		if ene_max_string[0] == '0' or  ene_max_string[0] == '.':
			if flag == 1:
				ene_max_string_notzero = ene_max_string_notzero[:5]
		else:
			flag = 0
		ene_max_string = ene_max_string_notzero	
	
	ene_type = str(ene_min_string)+'-'+str(ene_max_string)
	


if ene_range == 3:
	ene_dis = 'LIN'
	
	ene_min_string = repr(ene_min)	
	if type(ene_min) is int:
		pass
	else:
		nstring = len(ene_min_string)
		ene_min_string_notzero = ene_min_string
      		flag = 1 	
		if ene_min_string[0] == '0' or  ene_min_string[0] == '.':
			if flag == 1:
				ene_min_string_notzero = ene_min_string_notzero[:5]
		else:
			flag = 0
		ene_min_string = ene_min_string_notzero

	
	ene_max_string = repr(ene_max)	
	if type(ene_max) is int:
		pass
	else:
		nstring = len(ene_max_string)
		ene_max_string_notzero = ene_max_string
      		flag = 1 	
		if ene_max_string[0] == '0' or  ene_max_string[0] == '.':
			if flag == 1:
				ene_max_string_notzero = ene_max_string_notzero[:5]
		else:
			flag = 0
		ene_max_string = ene_max_string_notzero	
	
	ene_type = str(ene_min_string)+'-'+str(ene_max_string)

#print(ene_type)

if py_list == 0: 
	py_dir = 'QGSP_BERT_EMV'
	py_name = 'QGSP_BERT_EMV'

if py_list == 100:
	py_dir = '100List'
	py_name = '100List'

if py_list == 300:
	py_dir = '300List'
	py_name = '300List'

if py_list == 400:
	py_dir = 'ASTROMEV'
	py_name = 'ASTROMEV'


if sim_type == 0:
	sim_name = 'MONO'

if sim_type == 1:
	sim_name = 'RANGE'

if sim_type == 2:
	sim_name = 'CHEN'

if sim_type == 3:
	sim_name = 'VELA'

if sim_type == 4:
	sim_name = 'CRAB'

if sim_type == 5:
	sim_name = 'G400'


if pol_type == 1:
	pol_string = str(repr(pol_angle))+'POL.'
else: 
	pol_string = ''


if source_g == 0:
	sdir = '/Point'
	sname = 'Point'

if source_g == 1:
	sdir = '/Plane'
	sname = 'Plane'

#print(sdir, pol_string)


if cal_flag == 0 and ac_flag == 0:
	dir_cal = '/OnlyTracker'
if cal_flag == 1 and ac_flag == 0:
	dir_cal = '/onlyCAL'
if cal_flag == 0 and ac_flag == 1:
	dir_cal = '/onlyAC'
if cal_flag == 1 and ac_flag == 1:
	dir_cal = ''

if passive_flag == 0:
	dir_passive = ''
if passive_flag == 1:
	dir_passive = '/WithPassive'

if astrogam_version == 'V1.0':
	if isStrip == 0:
		stripDir = 'NoPixel/'
	if isStrip == 1 and repli == 0:
		stripDir = 'PixelNoRepli/'
	if isStrip == 1 and repli == 1:
		stripDir = 'PixelRepli/'
    
	if isStrip == 0:
		stripname = 'NOPIXEL'
	if isStrip == 1 and repli == 0:
		stripname = 'PIXEL'
	if isStrip == 1 and repli == 1:
		stripname = 'PIXEL.REPLI'


n_fits = os.listdir('./eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+sim_name+'/'+ene_type+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV/')

filepath = './eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+sim_name+'/'+ene_type+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV/'
print('LEVEL0 file path: '+ filepath)


aa_strip_event_id_list = []
aa_strip_theta_in_list = []
aa_strip_phi_in_list = []
aa_strip_ene_in_list = []
aa_strip_plane_id_list = []
aa_strip_zpos_list = []
aa_strip_si_id_list = []
aa_strip_strip_id_list = []
aa_strip_pos_list = []
aa_strip_edep_list = []
aa_strip_pair_list = []


aa_kalman_event_id_list = []
aa_kalman_theta_in_list = []
aa_kalman_phi_in_list = []
aa_kalman_ene_in_list = []
aa_kalman_plane_id_list = []
aa_kalman_zpos_list = []
aa_kalman_si_id_list = []
aa_kalman_pos_list = []
aa_kalman_edep_list = []
aa_kalman_strip_number_list = []
aa_kalman_pair_list = []


aa_kalman_pair_event_id_list = []
aa_kalman_pair_theta_in_list = []
aa_kalman_pair_phi_in_list = []
aa_kalman_pair_ene_in_list = []
aa_kalman_pair_plane_id_list = []
aa_kalman_pair_zpos_list = []
aa_kalman_pair_si_id_list = []
aa_kalman_pair_pos_list = []
aa_kalman_pair_edep_list = []
aa_kalman_pair_strip_number_list = []
aa_kalman_pair_pair_list = []


aa_kalman_compton_event_id_list = []
aa_kalman_compton_theta_in_list = []
aa_kalman_compton_phi_in_list = []
aa_kalman_compton_ene_in_list = []
aa_kalman_compton_plane_id_list = []
aa_kalman_compton_zpos_list = []
aa_kalman_compton_si_id_list = []
aa_kalman_compton_pos_list = []
aa_kalman_compton_edep_list = []
aa_kalman_compton_strip_number_list = []
aa_kalman_compton_pair_list = []


rawData_event_id_list = []
rawData_tray_id_list = []
rawData_plane_id_list = []
rawData_Strip_id_x_list = []
rawData_Strip_id_y_list = []
rawData_energy_dep_list = []
rawData_ent_x_list = []
rawData_ent_y_list = []
rawData_ent_z_list = []
rawData_exit_x_list = []
rawData_exit_y_list = []
rawData_exit_z_list = []
rawData_child_id_list = []
rawData_proc_id_list = []


L0TRACKER_Glob_event_id_list = []
L0TRACKER_Glob_vol_id_list = []		
L0TRACKER_Glob_moth_id_list = []		
L0TRACKER_Glob_tray_id_list = []
L0TRACKER_Glob_plane_id_list = []
L0TRACKER_Glob_Si_id_list = []
L0TRACKER_Glob_Strip_id_list = []
L0TRACKER_Glob_pos_list = []
L0TRACKER_Glob_zpos_list = []
L0TRACKER_Glob_energy_dep_list = []
L0TRACKER_Glob_pair_flag_list = []


L05TRACKER_Glob_event_id_cluster_list = []
L05TRACKER_Glob_tray_id_cluster_list = []
L05TRACKER_Glob_plane_id_cluster_list = []
L05TRACKER_Glob_Si_id_cluster_list = []
L05TRACKER_Glob_pos_cluster_list = []
L05TRACKER_Glob_zpos_cluster_list = []
L05TRACKER_Glob_energy_cluster_dep_list = []
L05TRACKER_Glob_pair_flag_cluster_list = []


aa_fake_event_id_list = []
aa_fake_theta_in_list = []
aa_fake_phi_in_list = []
aa_fake_ene_in_list = []
aa_fake_plane_id_list = []
aa_fake_zpos_list = []
aa_fake_si_id_list = []
aa_fake_pos_list = []
aa_fake_edep_list = []
aa_fake_strip_number_list = []
aa_fake_child_id_list = []
aa_fake_proc_id_list = []


calInput_event_id_tot_cal_list = []
calInput_bar_id_tot_list = []
calInput_bar_ene_tot_list = []


calInputSum_event_id_tot_cal_list = []
calInputSum_bar_ene_tot_list = []


acInput_event_id_tot_ac_list = []
acInput_AC_panel_list = []
acInput_AC_subpanel_list = []
acInput_energy_dep_tot_ac_list = []

############################################### Creazione file

ifile = 0
while ifile < len(n_fits):

	if isStrip == 1:
	
		filenamedat_aa_strip = filepath+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'

		table_aa_strip = open(filenamedat_aa_strip)
		
		lst = []
		
		for line in table_aa_strip:
			lst += [line.split()]
		
		event_ID = [x[0] for x in lst]
		theta_input = [x[1] for x in lst]
		phi_input = [x[2] for x in lst]
		energy_input = [x[3] for x in lst]
		plane_ID = [x[4] for x in lst]
		Pos_Z = [x[5] for x in lst]
		X_Y_flag = [x[6] for x in lst]
		strip_ID = [x[7] for x in lst]
		Strip_position = [x[8] for x in lst]
		energy_deposition = [x[9] for x in lst]
		pair_flag = [x[10] for x in lst]
	
		aa_strip_event_id_list.append(event_ID)
		aa_strip_theta_in_list.append(theta_input)
		aa_strip_phi_in_list.append(phi_input)
		aa_strip_ene_in_list.append(energy_input)
		aa_strip_plane_id_list.append(plane_ID)
		aa_strip_zpos_list.append(Pos_Z)
		aa_strip_si_id_list.append(X_Y_flag)
		aa_strip_strip_id_list.append(strip_ID)
		aa_strip_pos_list.append(Strip_position)
		aa_strip_edep_list.append(energy_deposition)
		aa_strip_pair_list.append(pair_flag)
		
		
		table_aa_strip.close()



		filenamedat_aa_kalman = filepath+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'
		
		table_aa_kalman_cluster = open(filenamedat_aa_kalman)
		
		lst = []
		
		for line in table_aa_strip:
			lst += [line.split()]
		
		event_ID = [x[0] for x in lst]
		theta_input = [x[1] for x in lst]
		phi_input = [x[2] for x in lst]
		energy_input = [x[3] for x in lst]
		plane_ID = [x[4] for x in lst]
		Pos_Z = [x[5] for x in lst]
		X_Y_flag = [x[6] for x in lst]
		Cluster_position = [x[7] for x in lst]
		energy_deposition = [x[8] for x in lst]
		number_of_strips_composing_the_cluster = [x[9] for x in lst]
		pair_flag = [x[10] for x in lst]
		
		aa_kalman_event_id_list.append(event_ID)
		aa_kalman_theta_in_list.append(theta_input)
		aa_kalman_phi_in_list.append(phi_input)
		aa_kalman_ene_in_list.append(energy_input)
		aa_kalman_plane_id_list.append(plane_ID)
		aa_kalman_zpos_list.append(Pos_Z)
		aa_kalman_si_id_list.append(X_Y_flag)
		aa_kalman_pos_list.append(Cluster_position)
		aa_kalman_edep_list.append(energy_deposition)
		aa_kalman_strip_number_list.append(number_of_strips_composing_the_cluster)
		aa_kalman_pair_list.append(pair_flag)
		
		
		table_aa_kalman_cluster.close()


		filenamedat_aa_kalman_pair = filepath+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'
		
		table_aa_kalman_cluster_pair = open(filenamedat_aa_kalman_pair)
		
		lst = []
		
		for line in table_aa_strip:
			lst += [line.split()]
		
		event_ID = [x[0] for x in lst]
		theta_input = [x[1] for x in lst]
		phi_input = [x[2] for x in lst]
		energy_input = [x[3] for x in lst]
		plane_ID = [x[4] for x in lst]
		Pos_Z = [x[5] for x in lst]
		X_Y_flag = [x[6] for x in lst]
		Cluster_position = [x[7] for x in lst]
		energy_deposition = [x[8] for x in lst]
		number_of_strips_composing_the_cluster = [x[9] for x in lst]
		pair_flag = [x[10] for x in lst]
		
		aa_kalman_pair_event_id_list.append(event_ID)
		aa_kalman_pair_theta_in_list.append(theta_input)
		aa_kalman_pair_phi_in_list.append(phi_input)
		aa_kalman_pair_ene_in_list.append(energy_input)
		aa_kalman_pair_plane_id_list.append(plane_ID)
		aa_kalman_pair_zpos_list.append(Pos_Z)
		aa_kalman_pair_si_id_list.append(X_Y_flag)
		aa_kalman_pair_pos_list.append(Cluster_position)
		aa_kalman_pair_edep_list.append(energy_deposition)
		aa_kalman_pair_strip_number_list.append(number_of_strips_composing_the_cluster)
		aa_kalman_pair_pair_list.append(pair_flag)
		
		
		table_aa_kalman_cluster_pair.close()


		filenamedat_aa_kalman_compton = filepath+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'
		
		table_aa_kalman_cluster_compton = open(filenamedat_aa_kalman_compton)
		
		lst = []
		
		for line in table_aa_strip:
			lst += [line.split()]
		
		event_ID = [x[0] for x in lst]
		theta_input = [x[1] for x in lst]
		phi_input = [x[2] for x in lst]
		energy_input = [x[3] for x in lst]
		plane_ID = [x[4] for x in lst]
		Pos_Z = [x[5] for x in lst]
		X_Y_flag = [x[6] for x in lst]
		Cluster_position = [x[7] for x in lst]
		energy_deposition = [x[8] for x in lst]
		number_of_strips_composing_the_cluster = [x[9] for x in lst]
		pair_flag = [x[10] for x in lst]
		
		aa_kalman_compton_event_id_list.append(event_ID)
		aa_kalman_compton_theta_in_list.append(theta_input)
		aa_kalman_compton_phi_in_list.append(phi_input)
		aa_kalman_compton_ene_in_list.append(energy_input)
		aa_kalman_compton_plane_id_list.append(plane_ID)
		aa_kalman_compton_zpos_list.append(Pos_Z)
		aa_kalman_compton_si_id_list.append(X_Y_flag)
		aa_kalman_compton_pos_list.append(Cluster_position)
		aa_kalman_compton_edep_list.append(energy_deposition)
		aa_kalman_compton_strip_number_list.append(number_of_strips_composing_the_cluster)
		aa_kalman_compton_pair_list.append(pair_flag)
		
		
		table_aa_kalman_cluster_compton.close()


		filenamefits_raw = fits.open(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_raw[1].data

		rawData_event_id_temp = tbdata.field('EVT_ID')
		rawData_tray_id_temp = tbdata.field('TRAY_ID')
		rawData_plane_id_temp = tbdata.field('PLANE_ID')
		rawData_Strip_id_x_temp = tbdata.field('STRIP_ID_X')
		rawData_Strip_id_y_temp = tbdata.field('STRIP_ID_Y')
		rawData_energy_dep_temp = tbdata.field('E_DEP')
		rawData_ent_x_temp = tbdata.field('X_ENT')
		rawData_ent_y_temp = tbdata.field('Y_ENT')
		rawData_ent_z_temp = tbdata.field('Z_ENT')
		rawData_exit_x_temp = tbdata.field('X_EXIT')
		rawData_exit_y_temp = tbdata.field('Y_EXIT')
		rawData_exit_z_temp = tbdata.field('Z_EXIT')
		rawData_child_id_temp = tbdata.field('CHILD_ID')
		rawData_proc_id_temp = tbdata.field('PROC_ID')


		rawData_event_id_list.append(rawData_event_id_temp)
		rawData_tray_id_list.append(rawData_tray_id_temp)
		rawData_plane_id_list.append(rawData_plane_id_temp)
		rawData_Strip_id_x_list.append(rawData_Strip_id_x_temp)
		rawData_Strip_id_y_list.append(rawData_Strip_id_y_temp)
		rawData_energy_dep_list.append(rawData_energy_dep_temp)
		rawData_ent_x_list.append(rawData_ent_x_temp)
		rawData_ent_y_list.append(rawData_ent_y_temp)
		rawData_ent_z_list.append(rawData_ent_z_temp)
		rawData_exit_x_list.append(rawData_exit_x_temp)
		rawData_exit_y_list.append(rawData_exit_y_temp)
		rawData_exit_z_list.append(rawData_exit_z_temp)
		rawData_child_id_list.append(rawData_child_id_temp)
		rawData_proc_id_list.append(rawData_proc_id_temp)

		filenamefits_raw.close()
		

		filenamefits_l0 = fits.open(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_l0[1].data

		L0TRACKER_Glob_event_id_temp = tbdata.field('EVT_ID')
		L0TRACKER_Glob_vol_id_temp = tbdata.field('VOLUME_ID')		
		L0TRACKER_Glob_moth_id_temp = tbdata.field('MOTHER_ID')		
		L0TRACKER_Glob_tray_id_temp = tbdata.field('TRAY_ID')
		L0TRACKER_Glob_plane_id_temp = tbdata.field('PLANE_ID')
		L0TRACKER_Glob_Si_id_temp = tbdata.field('TRK_FLAG')
		L0TRACKER_Glob_Strip_id_temp = tbdata.field('STRIP_ID')
		L0TRACKER_Glob_pos_temp = tbdata.field('POS')
		L0TRACKER_Glob_zpos_temp = tbdata.field('Z_POS')
		L0TRACKER_Glob_energy_dep_temp = tbdata.field('E_DEP')
		L0TRACKER_Glob_pair_flag_temp = tbdata.field('PAIR_FLAG')
		

		L0TRACKER_Glob_event_id_list.append(L0TRACKER_Glob_event_id_temp)
		L0TRACKER_Glob_vol_id_list.append(L0TRACKER_Glob_vol_id_temp)
		L0TRACKER_Glob_moth_id_list.append(L0TRACKER_Glob_moth_id_temp)		
		L0TRACKER_Glob_tray_id_list.append(L0TRACKER_Glob_tray_id_temp)
		L0TRACKER_Glob_plane_id_list.append(L0TRACKER_Glob_plane_id_temp)
		L0TRACKER_Glob_Si_id_list.append(L0TRACKER_Glob_Si_id_temp)
		L0TRACKER_Glob_Strip_id_list.append(L0TRACKER_Glob_Strip_id_temp)
		L0TRACKER_Glob_pos_list.append(L0TRACKER_Glob_pos_temp)
		L0TRACKER_Glob_zpos_list.append(L0TRACKER_Glob_zpos_temp)
		L0TRACKER_Glob_energy_dep_list.append(L0TRACKER_Glob_energy_dep_temp)
		L0TRACKER_Glob_pair_flag_list.append(L0TRACKER_Glob_pair_flag_temp)

		filenamefits_l0.close()


		filenamefits_l05 = fits.open(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_l05[1].data

		L05TRACKER_Glob_event_id_cluster_temp = tbdata.field('EVT_ID')
		L05TRACKER_Glob_tray_id_cluster_temp = tbdata.field('TRAY_ID')
		L05TRACKER_Glob_plane_id_cluster_temp = tbdata.field('PLANE_ID')
		L05TRACKER_Glob_Si_id_cluster_temp = tbdata.field('TRK_FLAG')
		L05TRACKER_Glob_pos_cluster_temp = tbdata.field('POS')
		L05TRACKER_Glob_zpos_cluster_temp = tbdata.field('Z_POS')
		L05TRACKER_Glob_energy_dep_cluster_temp = tbdata.field('E_DEP')
		L05TRACKER_Glob_pair_flag_cluster_temp = tbdata.field('PAIR_FLAG')
		

		L05TRACKER_Glob_event_id_cluster_list.append(L05TRACKER_Glob_event_id_cluster_temp)
		L05TRACKER_Glob_tray_id_cluster_list.append(L05TRACKER_Glob_tray_id_cluster_temp)
		L05TRACKER_Glob_plane_id_cluster_list.append(L05TRACKER_Glob_plane_id_cluster_temp)
		L05TRACKER_Glob_Si_id_cluster_list.append(L05TRACKER_Glob_Si_id_cluster_temp)
		L05TRACKER_Glob_pos_cluster_list.append(L05TRACKER_Glob_pos_cluster_temp)
		L05TRACKER_Glob_zpos_cluster_list.append(L05TRACKER_Glob_zpos_cluster_temp)
		L05TRACKER_Glob_energy_cluster_dep_list.append(L05TRACKER_Glob_energy_dep_cluster_temp)
		L05TRACKER_Glob_pair_flag_cluster_list.append(L05TRACKER_Glob_pair_flag_cluster_temp)

		filenamefits_l05.close()		

	else:

		filenamedat_aa_fake = filepath+'AA_FAKE_eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'
		
		table_aa_fake = ascii.open(filenamedat_aa_fake)
		
		lst = []
		
		for line in table_aa_strip:
			lst += [line.split()]
		
		event_ID = [x[0] for x in lst]
		theta_input = [x[1] for x in lst]
		phi_input = [x[2] for x in lst]
		energy_input = [x[3] for x in lst]
		plane_ID = [x[4] for x in lst]
		Pos_Z = [x[5] for x in lst]
		X_Y_flag = [x[6] for x in lst]
		Cluster_position = [x[7] for x in lst]
		energy_deposition = [x[8] for x in lst]
		number_of_strips_composing_the_cluster = [x[9] for x in lst]
		child_id = [x[10] for x in lst]
		proc_id = [x[11] for x in lst]
		
		aa_fake_event_id_list.append(event_ID)
		aa_fake_theta_in_list.append(theta_input)
		aa_fake_phi_in_list.append(phi_input)
		aa_fake_ene_in_list.append(energy_input)
		aa_fake_plane_id_list.append(plane_ID)
		aa_fake_zpos_list.append(Pos_Z)
		aa_fake_si_id_list.append(X_Y_flag)
		aa_fake_pos_list.append(Cluster_position)
		aa_fake_edep_list.append(energy_deposition)
		aa_fake_strip_number_list.append(number_of_strips_composing_the_cluster)
		aa_fake_child_id_list.append(child_id)
		aa_fake_proc_id_list.append(proc_id)
		
		
		table_aa_fake.close()


	if cal_flag == 1:
		
		filenamefits_cal = fits.open(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_cal[1].data

		calInput_event_id_tot_cal_temp = tbdata.field('EVT_ID')
		calInput_bar_id_tot_temp = tbdata.field('BAR_ID')
		calInput_bar_ene_tot_temp = tbdata.field('BAR_ENERGY')

		calInput_event_id_tot_cal_list.append(calInput_event_id_tot_cal_temp)
		calInput_bar_id_tot_list.append(calInput_bar_id_tot_temp)
		calInput_bar_ene_tot_list.append(calInput_bar_ene_tot_temp)

		filenamefits_cal.close()
	

		filenamefits_cal_sum = fits.open(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_cal_sum[1].data

		calInputSum_event_id_tot_cal_temp = tbdata.field('EVT_ID')
		calInputSum_bar_ene_tot_temp = tbdata.field('BAR_ENERGY')

		calInputSum_event_id_tot_cal_list.append(calInputSum_event_id_tot_cal_temp)
		calInputSum_bar_ene_tot_list.append(calInputSum_bar_ene_tot_temp)

		filenamefits_cal_sum.close()


	if ac_flag == 1:
		
		filenamefits_ac = fits.open(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_ac[1].data

		acInput_event_id_tot_ac_temp = tbdata.field('EVT_ID')
		acInput_AC_panel_temp = tbdata.field('AC_PANEL')
		acInput_AC_subpanel_temp = tbdata.field('AC_SUBPANEL')
		acInput_energy_dep_tot_ac_temp = tbdata.field('E_DEP')

		acInput_event_id_tot_ac_list.append(acInput_event_id_tot_ac_temp)
		acInput_AC_panel_list.append(acInput_AC_panel_temp)
		acInput_AC_subpanel_list.append(acInput_AC_subpanel_temp)
		acInput_energy_dep_tot_ac_list.append(acInput_energy_dep_tot_ac_temp)

		filenamefits_ac.close()		
#############################################

	ifile = ifile + 1

aa_strip_event_id = np.array(aa_strip_event_id_list)
aa_strip_theta_in = np.array(aa_strip_theta_in_list)
aa_strip_phi_in = np.array(aa_strip_phi_in_list)
aa_strip_ene_in = np.array(aa_strip_ene_in_list)
aa_strip_plane_id = np.array(aa_strip_plane_id_list)
aa_strip_zpos = np.array(aa_strip_zpos_list)
aa_strip_si_id = np.array(aa_strip_si_id_list)
aa_strip_strip_id = np.array(aa_strip_strip_id_list)
aa_strip_pos = np.array(aa_strip_pos_list)
aa_strip_edep = np.array(aa_strip_edep_list)
aa_strip_pair = np.array(aa_strip_pair_list)


aa_kalman_event_id = np.array(aa_kalman_event_id_list)
aa_kalman_theta_in = np.array(aa_kalman_theta_in_list)
aa_kalman_phi_in = np.array(aa_kalman_phi_in_list)
aa_kalman_ene_in = np.array(aa_kalman_ene_in_list)
aa_kalman_plane_id = np.array(aa_kalman_plane_id_list)
aa_kalman_zpos = np.array(aa_kalman_zpos_list)
aa_kalman_si_id = np.array(aa_kalman_si_id_list)
aa_kalman_pos = np.array(aa_kalman_pos_list)
aa_kalman_edep = np.array(aa_kalman_edep_list)
aa_kalman_strip_number = np.array(aa_kalman_strip_number_list)
aa_kalman_pair = np.array(aa_kalman_pair_list)


aa_kalman_pair_event_id = np.array(aa_kalman_pair_event_id_list)
aa_kalman_pair_theta_in = np.array(aa_kalman_pair_theta_in_list)
aa_kalman_pair_phi_in = np.array(aa_kalman_pair_phi_in_list)
aa_kalman_pair_ene_in = np.array(aa_kalman_pair_ene_in_list)
aa_kalman_pair_plane_id = np.array(aa_kalman_pair_plane_id_list)
aa_kalman_pair_zpos = np.array(aa_kalman_pair_zpos_list)
aa_kalman_pair_si_id = np.array(aa_kalman_pair_si_id_list)
aa_kalman_pair_pos = np.array(aa_kalman_pair_pos_list)
aa_kalman_pair_edep = np.array(aa_kalman_pair_edep_list)
aa_kalman_pair_strip_number = np.array(aa_kalman_pair_strip_number_list)
aa_kalman_pair_pair = np.array(aa_kalman_pair_pair_list)


aa_kalman_compton_event_id = np.array(aa_kalman_compton_event_id_list)
aa_kalman_compton_theta_in = np.array(aa_kalman_compton_theta_in_list)
aa_kalman_compton_phi_in = np.array(aa_kalman_compton_phi_in_list)
aa_kalman_compton_ene_in = np.array(aa_kalman_compton_ene_in_list)
aa_kalman_compton_plane_id = np.array(aa_kalman_compton_plane_id_list)
aa_kalman_compton_zpos = np.array(aa_kalman_compton_zpos_list)
aa_kalman_compton_si_id = np.array(aa_kalman_compton_si_id_list)
aa_kalman_compton_pos = np.array(aa_kalman_compton_pos_list)
aa_kalman_compton_edep = np.array(aa_kalman_compton_edep_list)
aa_kalman_compton_strip_number = np.array(aa_kalman_compton_strip_number_list)
aa_kalman_compton_pair = np.array(aa_kalman_compton_pair_list)


rawData_event_id = np.array(rawData_event_id_list)
rawData_tray_id = np.array(rawData_tray_id_list)
rawData_plane_id = np.array(rawData_plane_id_list)
rawData_Strip_id_x = np.array(rawData_Strip_id_x_list)
rawData_Strip_id_y = np.array(rawData_Strip_id_y_list)
rawData_energy_dep = np.array(rawData_energy_dep_list)
rawData_ent_x = np.array(rawData_ent_x_list)
rawData_ent_y = np.array(rawData_ent_y_list)
rawData_ent_z = np.array(rawData_ent_z_list)
rawData_exit_x = np.array(rawData_exit_x_list)
rawData_exit_y = np.array(rawData_exit_y_list)
rawData_exit_z = np.array(rawData_exit_z_list)
rawData_child_id = np.array(rawData_child_id_list)
rawData_proc_id = np.array(rawData_proc_id_list)


L0TRACKER_Glob_event_id = np.array(L0TRACKER_Glob_event_id_list)
L0TRACKER_Glob_vol_id = np.array(L0TRACKER_Glob_vol_id_list)		
L0TRACKER_Glob_moth_id = np.array(L0TRACKER_Glob_moth_id_list)		
L0TRACKER_Glob_tray_id = np.array(L0TRACKER_Glob_tray_id_list)
L0TRACKER_Glob_plane_id = np.array(L0TRACKER_Glob_plane_id_list)
L0TRACKER_Glob_Si_id = np.array(L0TRACKER_Glob_Si_id_list)
L0TRACKER_Glob_Strip_id = np.array(L0TRACKER_Glob_Strip_id_list)
L0TRACKER_Glob_pos = np.array(L0TRACKER_Glob_pos_list)
L0TRACKER_Glob_zpos = np.array(L0TRACKER_Glob_zpos_list)
L0TRACKER_Glob_energy_dep = np.array(L0TRACKER_Glob_energy_dep_list)
L0TRACKER_Glob_pair_flag = np.array(L0TRACKER_Glob_pair_flag_list)


L05TRACKER_Glob_event_id_cluster = np.array(L05TRACKER_Glob_event_id_cluster_list)
L05TRACKER_Glob_tray_id_cluster = np.array(L05TRACKER_Glob_tray_id_cluster_list)
L05TRACKER_Glob_plane_id_cluster = np.array(L05TRACKER_Glob_plane_id_cluster_list)
L05TRACKER_Glob_Si_id_cluster = np.array(L05TRACKER_Glob_Si_id_cluster_list)
L05TRACKER_Glob_pos_cluster = np.array(L05TRACKER_Glob_pos_cluster_list)
L05TRACKER_Glob_zpos_cluster = np.array(L05TRACKER_Glob_zpos_cluster_list)
L05TRACKER_Glob_energy_cluster_dep = np.array(L05TRACKER_Glob_energy_cluster_dep_list)
L05TRACKER_Glob_pair_flag_cluster = np.array(L05TRACKER_Glob_pair_flag_cluster_list)


aa_fake_event_id = np.array(aa_fake_event_id_list)
aa_fake_theta_in = np.array(aa_fake_theta_in_list)
aa_fake_phi_in = np.array(aa_fake_phi_in_list)
aa_fake_ene_in = np.array(aa_fake_ene_in_list)
aa_fake_plane_id = np.array(aa_fake_plane_id_list)
aa_fake_zpos = np.array(aa_fake_zpos_list)
aa_fake_si_id = np.array(aa_fake_si_id_list)
aa_fake_pos = np.array(aa_fake_pos_list)
aa_fake_edep = np.array(aa_fake_edep_list)
aa_fake_strip_number = np.array(aa_fake_strip_number_list)
aa_fake_child_id = np.array(aa_fake_child_id_list)
aa_fake_proc_id = np.array(aa_fake_proc_id_list)


calInput_event_id_tot_cal = np.array(calInput_event_id_tot_cal_list)
calInput_bar_id_tot = np.array(calInput_bar_id_tot_list)
calInput_bar_ene_tot = np.array(calInput_bar_ene_tot_list)


calInputSum_event_id_tot_cal = np.array(calInputSum_event_id_tot_cal_list)
calInputSum_bar_ene_tot = np.array(calInputSum_bar_ene_tot_list)


acInput_event_id_tot_ac = np.array(acInput_event_id_tot_ac_list)
acInput_AC_panel = np.array(acInput_AC_panel_list)
acInput_AC_subpanel = np.array(acInput_AC_subpanel_list)
acInput_energy_dep_tot_ac = np.array(acInput_energy_dep_tot_ac_list)

############################### Scrittura file


if isStrip == 1:

	if os.path.exists(filepath+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat'):
		os.remove(filepath+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat')
		data = open(filepath+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')
	else:
		data = open(filepath+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')

	for r in range(len(aa_strip_event_id)):
		
		data.write('{:d}\t'.format(aa_strip_event_id[r]))
		data.write('{:d}\t'.format(aa_strip_theta_in[r]))
		data.write('{:d}\t'.format(aa_strip_phi_in[r]))
		data.write('{:s}\t'.format(aa_strip_ene_in[r]))
		data.write('{:d}\t'.format(aa_strip_plane_id[r]))
		data.write('{:f}\t'.format(aa_strip_zpos[r]))
		data.write('{:d}\t'.format(aa_strip_si_id[r]))
		data.write('{:f}\t'.format(aa_strip_strip_id[r]))
		data.write('{:f}\t'.format(aa_strip_pos[r]))
		data.write('{:d}\t'.format(aa_strip_edep[r]))
		data.write('{:d}\n'.format(aa_strip_pair[r]))

	data.close()
	
	
	if os.path.exists(filepath+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat'):
		os.remove(filepath+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat')
		data = open(filepath+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')
	else:
		data = open(filepath+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')

	for r in range(len(aa_kalman_event_id)):
		
		data.write('{:d}\t'.format(aa_kalman_event_id[r]))
		data.write('{:d}\t'.format(aa_kalman_theta_in[r]))
		data.write('{:d}\t'.format(aa_kalman_phi_in[r]))
		data.write('{:s}\t'.format(aa_kalman_ene_in[r]))
		data.write('{:d}\t'.format(aa_kalman_plane_id[r]))
		data.write('{:f}\t'.format(aa_kalman_zpos[r]))
		data.write('{:d}\t'.format(aa_kalman_si_id[r]))
		data.write('{:f}\t'.format(aa_kalman_pos[r]))
		data.write('{:f}\t'.format(aa_kalman_edep[r]))
		data.write('{:d}\t'.format(aa_kalman_strip_number[r]))
		data.write('{:d}\n'.format(aa_kalman_pair[r]))

	data.close()
	

	if os.path.exists(filepath+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat'):
		os.remove(filepath+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat')
		data = open(filepath+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')
	else:
		data = open(filepath+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')

	for r in range(len(aa_kalman_pair_event_id)):
		
		data.write('{:d}\t'.format(aa_kalman_pair_event_id[r]))
		data.write('{:d}\t'.format(aa_kalman_pair_theta_in[r]))
		data.write('{:d}\t'.format(aa_kalman_pair_phi_in[r]))
		data.write('{:s}\t'.format(aa_kalman_pair_ene_in[r]))
		data.write('{:d}\t'.format(aa_kalman_pair_plane_id[r]))
		data.write('{:f}\t'.format(aa_kalman_pair_zpos[r]))
		data.write('{:d}\t'.format(aa_kalman_pair_si_id[r]))
		data.write('{:f}\t'.format(aa_kalman_pair_pos[r]))
		data.write('{:f}\t'.format(aa_kalman_pair_edep[r]))
		data.write('{:d}\t'.format(aa_kalman_pair_strip_number[r]))
		data.write('{:d}\n'.format(aa_kalman_pair_pair[r]))

	data.close()


	if os.path.exists(filepath+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat'):
		os.remove(filepath+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat')
		data = open(filepath+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')
	else:
		data = open(filepath+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')

	for r in range(len(aa_kalman_compton_event_id)):
		
		data.write('{:d}\t'.format(aa_kalman_compton_event_id[r]))
		data.write('{:d}\t'.format(aa_kalman_compton_theta_in[r]))
		data.write('{:d}\t'.format(aa_kalman_compton_phi_in[r]))
		data.write('{:s}\t'.format(aa_kalman_compton_ene_in[r]))
		data.write('{:d}\t'.format(aa_kalman_compton_plane_id[r]))
		data.write('{:f}\t'.format(aa_kalman_compton_zpos[r]))
		data.write('{:d}\t'.format(aa_kalman_compton_si_id[r]))
		data.write('{:f}\t'.format(aa_kalman_compton_pos[r]))
		data.write('{:f}\t'.format(aa_kalman_compton_edep[r]))
		data.write('{:d}\t'.format(aa_kalman_compton_strip_number[r]))
		data.write('{:d}\n'.format(aa_kalman_compton_pair[r]))

	data.close()

else:

	if os.path.exists(filepath+'AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'_'+pol_string+'all.dat'):
		os.remove(filepath+'AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'_'+pol_string+'all.dat')
		data = open(filepath+'AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'_'+pol_string+'all.dat', 'w')
	else:
		data = open(filepath+'AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'_'+pol_string+'all.dat', 'w')

	for r in range(len(aa_fake_event_id)):
		
		data.write('{:d}\t'.format(aa_fake_event_id[r]))
		data.write('{:d}\t'.format(aa_fake_theta_in[r]))
		data.write('{:d}\t'.format(aa_fake_phi_in[r]))
		data.write('{:s}\t'.format(aa_fake_ene_in[r]))
		data.write('{:d}\t'.format(aa_fake_plane_id[r]))
		data.write('{:f}\t'.format(aa_fake_zpos[r]))
		data.write('{:d}\t'.format(aa_fake_si_id[r]))
		data.write('{:f}\t'.format(aa_fake_pos[r]))
		data.write('{:f}\t'.format(aa_fake_edep[r]))
		data.write('{:d}\t'.format(aa_fake_strip_number[r]))
		data.write('{:d}\t'.format(aa_fake_child_id[r]))
		data.write('{:d}\n'.format(aa_fake_proc_id[r]))

	data.close()


if isStrip == 1:

	col1 = fits.Column(name='EVT_ID', format='I', array=rawData_event_id)	
	col2 = fits.Column(name='TRAY_ID', format='I', array=rawData_tray_id)
	col3 = fits.Column(name='PLANE_ID', format='I', array=rawData_plane_id)
	col4 = fits.Column(name='STRIP_ID_X', format='I', array=rawData_Strip_id_x)
	col5 = fits.Column(name='STRIP_ID_Y', format='I', array=rawData_Strip_id_y)
	col6 = fits.Column(name='E_DEP', format='F20.5', array=rawData_energy_dep)
	col7 = fits.Column(name='X_ENT', format='F20.5', array=rawData_ent_x)
	col8 = fits.Column(name='Y_ENT', format='F20.5', array=rawData_ent_y)
	col9 = fits.Column(name='Z_ENT', format='F20.5', array=rawData_ent_z)
	col10 = fits.Column(name='X_EXIT', format='F20.5', array=rawData_exit_x)
	col11 = fits.Column(name='Y_EXIT', format='F20.5', array=rawData_exit_y)
	col12 = fits.Column(name='Z_EXIT', format='F20.5', array=rawData_exit_z)
	col13 = fits.Column(name='CHILD_ID', format='I', array=rawData_child_id)
	col14 = fits.Column(name='PROC_ID', format='I', array=rawData_proc_id)
		
	cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='CREATOR  = Giovanni Giannella & Simone Guidotti', ext=1)
		
	fits.setval(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='THELSIM release  = eASTROGAM '+astrogam_version, ext=1)
		
	fits.setval(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

	fits.setval(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='ENERGY ='+ene_type, ext=1)

	fits.setval(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)
		
	fits.setval(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

	fits.setval(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Position unit = cm', ext=1)
		
	fits.setval(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = keV', ext=1)




	col1 = fits.Column(name='EVT_ID', format='J', array=L0TRACKER_Glob_event_id)	
	col2 = fits.Column(name='VOLUME_ID', format='J', array=L0TRACKER_Glob_vol_id)
	col3 = fits.Column(name='MOTHER_ID', format='J', array=L0TRACKER_Glob_moth_id)
	col4 = fits.Column(name='TRAY_ID', format='I', array=L0TRACKER_Glob_tray_id)
	col5 = fits.Column(name='PLANE_ID', format='I', array=L0TRACKER_Glob_plane_id)
	col6 = fits.Column(name='TRK_FLAG', format='I', array=L0TRACKER_Glob_si_id)
	col7 = fits.Column(name='STRIP_ID', format='J', array=L0TRACKER_Glob_Strip_id)
	col8 = fits.Column(name='POS', format='F20.5', array=L0TRACKER_Glob_pos)
	col9 = fits.Column(name='ZPOS', format='F20.5', array=L0TRACKER_Glob_zpos)
	col10 = fits.Column(name='E_DEP', format='F20.5', array=L0TRACKER_Glob_energy_dep)
	col11 = fits.Column(name='PAIR_FLAG', format='I', array=L0TRACKER_Glob_pair_flag)

		
	cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='CREATOR  = Giovanni Giannella & Simone Guidotti', ext=1)

	fits.setval(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='THELSIM release  = eASTROGAM '+astrogam_version, ext=1)
		
	fits.setval(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in ='+str(N_in)+'      /Number of simulated particles', ext=1)

	fits.setval(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='ENERGY ='+ene_type+'     /Simulated input energy', ext=1)

	fits.setval(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta ='+str(theta_type)+'     /Simulated input theta angle', ext=1)
		
	fits.setval(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi ='+str(phi_type)+'     /Simulated input phi angle', ext=1)

	fits.setval(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = keV', ext=1)



	col1 = fits.Column(name='EVT_ID', format='J', array=L05TRACKER_Glob_event_id_cluster)	
	col2 = fits.Column(name='TRAY_ID', format='I', array=L05TRACKER_Glob_tray_id_cluster)
	col3 = fits.Column(name='PLANE_ID', format='I', array=L05TRACKER_Glob_plane_id_cluster)
	col4 = fits.Column(name='TRK_FLAG', format='I', array=L05TRACKER_Glob_si_id_cluster)
	col5 = fits.Column(name='POS', format='F20.5', array=L0TRACKER_Glob_pos_cluster)
	col6 = fits.Column(name='ZPOS', format='F20.5', array=L0TRACKER_Glob_zpos_cluster)
	col7 = fits.Column(name='E_DEP', format='F20.5', array=L0TRACKER_Glob_energy_dep_cluster)
	col8 = fits.Column(name='PAIR_FLAG', format='I', array=L0TRACKER_Glob_pair_flag_cluster)

		
	cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='CREATOR  = Giovanni Giannella & Simone Guidotti', ext=1)

	fits.setval(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='THELSIM release  = eASTROGAM '+astrogam_version, ext=1)
		
	fits.setval(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in ='+str(N_in)+'      /Number of simulated particles', ext=1)

	fits.setval(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='ENERGY ='+ene_type+'     /Simulated input energy', ext=1)

	fits.setval(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta ='+str(theta_type)+'     /Simulated input theta angle', ext=1)
		
	fits.setval(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi ='+str(phi_type)+'     /Simulated input phi angle', ext=1)

	fits.setval(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = keV', ext=1)


if cal_flag == 1:

	col1 = fits.Column(name='EVT_ID', format='I', array=calInput_event_id_tot_cal)	
	col2 = fits.Column(name='BAR_ID', format='I', array=calInput_bar_id_tot)
	col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=calInput_bar_ene_tot)
		
	cols = fits.ColDefs([col1,col2,col3])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='CREATOR  = Giovanni Giannella & Simone Guidotti', ext=1)

	fits.setval(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)
		
	fits.setval(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

	fits.setval(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='ENERGY ='+ene_type, ext=1)

	fits.setval(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)
		
	fits.setval(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

	fits.setval(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)
	

	col1 = fits.Column(name='EVT_ID', format='I', array=calInputSum_event_id_tot_cal)	
	col2 = fits.Column(name='BAR_ENERGY', format='F20.15', array=calInputSum_bar_ene_tot)
		
	cols = fits.ColDefs([col1,col2])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='CREATOR  = Giovanni Giannella & Simone Guidotti', ext=1)

	fits.setval(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)
		
	fits.setval(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

	fits.setval(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='ENERGY ='+ene_type, ext=1)

	fits.setval(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)
		
	fits.setval(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

	fits.setval(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)



if ac_flag == 1:

	col1 = fits.Column(name='EVT_ID', format='I', array=acInput_event_id_tot_ac)	
	col2 = fits.Column(name='AC_PANEL', format='A', array=acInput_AC_panel)
	col3 = fits.Column(name='AC_SUBPANEL', format='I', array=acInput_AC_subpanel)
	col4 = fits.Column(name='E_DEP', format='F20.15', array=acInput_energy_dep_tot_ac)
		
	cols = fits.ColDefs([col1,col2,col3,col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')


	fits.setval(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

	fits.setval(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

	fits.setval(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

	fits.setval(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

	fits.setval(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)
	















