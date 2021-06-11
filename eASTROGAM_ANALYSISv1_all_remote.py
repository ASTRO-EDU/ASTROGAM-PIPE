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
import sys

sys.argv[0] = 'eASTROGAM_ANALYSISv1_all_remote.py'

##############################
#     parametri iniziali     #
##############################

astrogam_version = sys.argv[1]           # Enter eASTROGAM release (e.g. V1.0, V2.0):
bogemms_tag = int(sys.argv[2])           # Enter BoGEMMS release (e.g. 211):
sim_type = int(sys.argv[3])              # Enter simulation type [0 = Mono, 1 = Range, 2 = Chen, 3: Vela, 4: Crab, 5: G400]:
py_list = int(sys.argv[4])               # Enter the Physics List [0 = QGSP_BERT_EMV, 100 = ARGO, 300 = FERMI, 400 = ASTROMEV]:
N_in = int(sys.argv[5])                  # Enter the number of emitted particles:
part_type = sys.argv[6]                  # Enter the particle type [ph = photons, mu = muons, g = geantino, p = proton, el = electron]:
ene_range = int(sys.argv[7])             # Enter energy distribution [0 = MONO, 1 = POW, 2 = EXP, 3 = LIN]:
ene_min = sys.argv[8]               	 # Enter miminum energy [MeV]:
ene_max = sys.argv[9]                    # Enter maximum energy [MeV]:
ang_type = sys.argv[10]                  # Enter the angular distribution [e.g. UNI, ISO]:
theta_type = int(sys.argv[11])           # Enter theta:
phi_type = int(sys.argv[12])             # Enter phi:
pol_type = int(sys.argv[13])             # Is the source polarized? [0 = false, 1 = true]:
pol_angle = int(sys.argv[14])            # Enter Polarization angle:
source_g = int(sys.argv[15])             # Enter source geometry [0 = Point, 1 = Plane]:
isStrip = int(sys.argv[16])              # Strip/Pixels activated? [0 = false, 1 = true]
repli = int(sys.argv[17])                # Strips/Pixels replicated? [0 = false, 1 = true]
cal_flag = int(sys.argv[18])             # Is Cal present? [0 = false, 1 = true]:
ac_flag = int(sys.argv[19])              # Is AC present? [0 = false, 1 = true]:
passive_flag = int(sys.argv[20])         # Is Passive present? [0 = false, 1 = true]:
energy_thresh = int(sys.argv[21])        # Enter energy threshold [keV]:
ifile = int(sys.argv[22])		 	  	 # Enter the initial number of FITS files:
n_fits = int(sys.argv[23])               # Enter the final number of FITS files:



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
	
if astrogam_version=='V1.1':
	astrogam_tag = '11'
	sim_tag = 'eAST'+str(bogemms_tag)+str(astrogam_tag)+'2021'

if astrogam_version=='V2.0':
    astrogam_tag = '20'
    sim_tag = 'eAST'+str(bogemms_tag)+str(astrogam_tag)+'2021'

if astrogam_version=='V10.0':
    astrogam_tag = '10'
    sim_tag = 'eAST'+str(bogemms_tag)+str(astrogam_tag)+'021'

if (ene_min[0] == '0'):
	ene_min = np.round(float(ene_min), 1)
else:
    ene_min_float = float(ene_min)
    if ene_min_float.is_integer():
        ene_min = int(ene_min_float)
    else:
        ene_min = np.round(ene_min_float, 1)

if (ene_max[0] == '0'):
	ene_max = np.round(float(ene_max), 1)
else:
    ene_max_float = float(ene_max)
    if ene_max_float.is_integer():
        ene_max = int(ene_max_float)
    else:
        ene_max = np.round(ene_max_float, 1)

if ene_range == 0:
	ene_dis = 'MONO'
	ene_type = str(ene_min)


if ene_range == 1:
	ene_dis = 'POW'

	ene_min_string = str(ene_min)	
	ene_max_string = str(ene_max)	

	ene_type = str(ene_min_string)+'.'+str(ene_max_string)




if ene_range == 2:
	ene_dis = 'EXP'

	ene_min_string = str(ene_min)	
	ene_max_string = str(ene_max)	

	ene_type = str(ene_min_string)+'.'+str(ene_max_string)



if ene_range == 3:
	ene_dis = 'LIN'

	ene_min_string = str(ene_min)	
	ene_max_string = str(ene_max)

	ene_type = str(ene_min_string)+'.'+str(ene_max_string)



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

if astrogam_version == 'V1.0' or astrogam_version == 'V1.1' or astrogam_version == 'V2.0' or astrogam_version == 'V10.0':
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

# READING THE FITS FILES

filepath = './output_eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+sim_name+'/'+ene_type+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV/'

outdir = ('./output_eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+str(sim_name)+'/'+str(ene_type)+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV')


print('LEVEL0 file path: '+ filepath)


# G4.RAW.TRACKER.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
rawData_event_id = []
rawData_vol_id = []
rawData_moth_id = []
rawData_tray_id = []
rawData_plane_id = []
rawData_Strip_id_x = []
rawData_Strip_id_y = []
rawData_energy_dep = []
rawData_ent_x = []
rawData_ent_y = []
rawData_ent_z = []
rawData_exit_x = []
rawData_exit_y = []
rawData_exit_z = []
rawData_part_id = []
rawData_trk_id = []
rawData_child_id = []
rawData_proc_id = []

if isStrip == 1:
	# L0.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	L0TRACKER_Glob_event_id = []
	L0TRACKER_Glob_vol_id = []		
	L0TRACKER_Glob_moth_id = []		
	L0TRACKER_Glob_tray_id = []
	L0TRACKER_Glob_plane_id = []
	L0TRACKER_Glob_Si_id = []
	L0TRACKER_Glob_Strip_id = []
	L0TRACKER_Glob_pos = []
	L0TRACKER_Glob_zpos = []
	L0TRACKER_Glob_energy_dep = []
	L0TRACKER_Glob_pair_flag = []

	# L0.5.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	L05TRACKER_Glob_event_id_cluster = []
	L05TRACKER_Glob_tray_id_cluster = []
	L05TRACKER_Glob_plane_id_cluster = []
	L05TRACKER_Glob_Si_id_cluster = []
	L05TRACKER_Glob_pos_cluster = []
	L05TRACKER_Glob_zpos_cluster = []
	L05TRACKER_Glob_energy_cluster_dep = []
	L05TRACKER_Glob_pair_flag_cluster = []

	# sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
	aa_strip_event_id = []
	aa_strip_theta_in = []
	aa_strip_phi_in = []
	aa_strip_ene_in = []
	aa_strip_plane_id = []
	aa_strip_zpos = []
	aa_strip_si_id = []
	aa_strip_strip_id = []
	aa_strip_pos = []
	aa_strip_edep = []
	aa_strip_pair = []

	# sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
	aa_kalman_event_id = []
	aa_kalman_theta_in = []
	aa_kalman_phi_in = []
	aa_kalman_ene_in = []
	aa_kalman_plane_id = []
	aa_kalman_zpos = []
	aa_kalman_si_id = []
	aa_kalman_pos = []
	aa_kalman_edep = []
	aa_kalman_strip_number = []
	aa_kalman_pair = []

	# sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
	aa_kalman_pair_event_id = []
	aa_kalman_pair_theta_in = []
	aa_kalman_pair_phi_in = []
	aa_kalman_pair_ene_in = []
	aa_kalman_pair_plane_id = []
	aa_kalman_pair_zpos = []
	aa_kalman_pair_si_id = []
	aa_kalman_pair_pos = []
	aa_kalman_pair_edep = []
	aa_kalman_pair_strip_number = []
	aa_kalman_pair_pair = []

	# sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
	aa_kalman_compton_event_id = []
	aa_kalman_compton_theta_in = []
	aa_kalman_compton_phi_in = []
	aa_kalman_compton_ene_in = []
	aa_kalman_compton_plane_id = []
	aa_kalman_compton_zpos = []
	aa_kalman_compton_si_id = []
	aa_kalman_compton_pos = []
	aa_kalman_compton_edep = []
	aa_kalman_compton_strip_number = []
	aa_kalman_compton_pair = []

	# sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
	aa_kalman_ray_event_id = []
	aa_kalman_ray_theta_in = []
	aa_kalman_ray_phi_in = []
	aa_kalman_ray_ene_in = []
	aa_kalman_ray_plane_id = []
	aa_kalman_ray_zpos = []
	aa_kalman_ray_si_id = []
	aa_kalman_ray_pos = []
	aa_kalman_ray_edep = []
	aa_kalman_ray_strip_number = []
	aa_kalman_ray_pair = []
	
	# sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'_TRACKER.dat')
	S1_event_id_trk = []
	S1_type_trk = []
	S1_e_dep_trk = []
	S1_volume_id_trk = []
	S1_x_trk = []
	S1_y_trk = []
	S1_z_trk = []
	
else:
	# sim_tag+'_AA_FAKE_eASTROGAM_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
	aa_fake_event_id = []
	aa_fake_theta_in = []
	aa_fake_phi_in = []
	aa_fake_ene_in = []
	aa_fake_plane_id = []
	aa_fake_zpos = []
	aa_fake_si_id = []
	aa_fake_pos = []
	aa_fake_edep = []
	aa_fake_strip_number = []
	aa_fake_child_id = []
	aa_fake_proc_id = []


if cal_flag == 1:
	# G4.RAW.CAL.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	rawData_event_id_cal = []
	rawData_vol_id_cal = []
	rawData_moth_id_cal = []
	rawData_energy_dep_cal = []
	rawData_ent_x_cal = []
	rawData_ent_y_cal = []
	rawData_ent_z_cal = []
	rawData_exit_x_cal = []
	rawData_exit_y_cal = []
	rawData_exit_z_cal = []
	rawData_part_id_cal = []
	rawData_trk_id_cal = []
	rawData_child_id_cal = []
	rawData_proc_id_cal = []

	# G4.CAL.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	calInput_event_id_tot = []
	calInput_bar_id_tot = []
	calInput_bar_ene_tot = []
	calInput_pair_flag_tot = []

	# G4.CAL.COMPTON.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	calInput_event_id_tot_compton = []
	calInput_bar_id_tot_compton = []
	calInput_bar_ene_tot_compton = []
	calInput_pair_flag_tot_compton = []

	# G4.CAL.PAIR.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	calInput_event_id_tot_pair = []
	calInput_bar_id_tot_pair = []
	calInput_bar_ene_tot_pair = []
	calInput_pair_flag_tot_pair = []

	# G4.CAL.RAYLEIGH.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	calInput_event_id_tot_ray = []
	calInput_bar_id_tot_ray = []
	calInput_bar_ene_tot_ray = []
	calInput_pair_flag_tot_ray = []
	
	# SUM.CAL.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	calInputSum_event_id_tot = []
	calInputSum_bar_ene_tot = []
	
	# sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'_CAL.dat')
	S1_event_id_cal = []
	S1_type_cal = []
	S1_e_dep_cal = []
	S1_volume_id_cal = []
	S1_x_cal = []
	S1_y_cal = []
	S1_z_cal = []


if ac_flag == 1:
	# G4.RAW.AC.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	rawData_event_id_ac = []
	rawData_vol_id_ac = []
	rawData_moth_id_ac = []
	rawData_energy_dep_ac = []
	rawData_ent_x_ac = []
	rawData_ent_y_ac = []
	rawData_ent_z_ac = []
	rawData_exit_x_ac = []
	rawData_exit_y_ac = []
	rawData_exit_z_ac = []
	rawData_part_id_ac = []
	rawData_trk_id_ac = []
	rawData_child_id_ac = []
	rawData_proc_id_ac = []

	# G4.AC.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	acInput_event_id_tot = []
	acInput_AC_panel = []
	acInput_AC_subpanel = []
	acInput_energy_dep_tot = []
	acInput_pair_flag_tot = []
	
	# G4.AC.COMPTON.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	acInput_event_id_tot_compton = []
	acInput_AC_panel_compton = []
	acInput_AC_subpanel_compton = []
	acInput_energy_dep_tot_compton = []
	acInput_pair_flag_tot_compton = []

	# G4.AC.PAIR.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	acInput_event_id_tot_pair = []
	acInput_AC_panel_pair = []
	acInput_AC_subpanel_pair = []
	acInput_energy_dep_tot_pair = []
	acInput_pair_flag_tot_pair = []

	# G4.AC.RAYLEIGH.eASTROGAM<version>.<phys>List.<strip>.<point>.<n_in>ph.<energy>MeV.<theta>.<phi>.all.fits
	acInput_event_id_tot_ray = []
	acInput_AC_panel_ray = []
	acInput_AC_subpanel_ray = []
	acInput_energy_dep_tot_ray = []
	acInput_pair_flag_tot_ray = []	
	
	# sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'_AC.dat')
	S1_event_id_ac = []
	S1_type_ac = []
	S1_e_dep_ac = []
	S1_volume_id_ac = []
	S1_x_ac = []
	S1_y_ac = []
	S1_z_ac = []	

############################################### Creazione file


while ifile <= n_fits:

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
	
		aa_strip_event_id.append(event_ID)
		aa_strip_theta_in.append(theta_input)
		aa_strip_phi_in.append(phi_input)
		aa_strip_ene_in.append(energy_input)
		aa_strip_plane_id.append(plane_ID)
		aa_strip_zpos.append(Pos_Z)
		aa_strip_si_id.append(X_Y_flag)
		aa_strip_strip_id.append(strip_ID)
		aa_strip_pos.append(Strip_position)
		aa_strip_edep.append(energy_deposition)
		aa_strip_pair.append(pair_flag)
		
		
		table_aa_strip.close()



		filenamedat_aa_kalman = filepath+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'
		
		table_aa_kalman_cluster = open(filenamedat_aa_kalman)
		
		lst = []
		
		for line in table_aa_kalman_cluster:
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
		
		aa_kalman_event_id.append(event_ID)
		aa_kalman_theta_in.append(theta_input)
		aa_kalman_phi_in.append(phi_input)
		aa_kalman_ene_in.append(energy_input)
		aa_kalman_plane_id.append(plane_ID)
		aa_kalman_zpos.append(Pos_Z)
		aa_kalman_si_id.append(X_Y_flag)
		aa_kalman_pos.append(Cluster_position)
		aa_kalman_edep.append(energy_deposition)
		aa_kalman_strip_number.append(number_of_strips_composing_the_cluster)
		aa_kalman_pair.append(pair_flag)
		
		
		table_aa_kalman_cluster.close()


		filenamedat_aa_kalman_pair = filepath+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'
		
		table_aa_kalman_cluster_pair = open(filenamedat_aa_kalman_pair)
		
		lst = []
		
		for line in table_aa_kalman_cluster_pair:
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
		
		aa_kalman_pair_event_id.append(event_ID)
		aa_kalman_pair_theta_in.append(theta_input)
		aa_kalman_pair_phi_in.append(phi_input)
		aa_kalman_pair_ene_in.append(energy_input)
		aa_kalman_pair_plane_id.append(plane_ID)
		aa_kalman_pair_zpos.append(Pos_Z)
		aa_kalman_pair_si_id.append(X_Y_flag)
		aa_kalman_pair_pos.append(Cluster_position)
		aa_kalman_pair_edep.append(energy_deposition)
		aa_kalman_pair_strip_number.append(number_of_strips_composing_the_cluster)
		aa_kalman_pair_pair.append(pair_flag)
		
		
		table_aa_kalman_cluster_pair.close()


		filenamedat_aa_kalman_compton = filepath+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'
		
		table_aa_kalman_cluster_compton = open(filenamedat_aa_kalman_compton)
		
		lst = []
		
		for line in table_aa_kalman_cluster_compton:
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
		
		aa_kalman_compton_event_id.append(event_ID)
		aa_kalman_compton_theta_in.append(theta_input)
		aa_kalman_compton_phi_in.append(phi_input)
		aa_kalman_compton_ene_in.append(energy_input)
		aa_kalman_compton_plane_id.append(plane_ID)
		aa_kalman_compton_zpos.append(Pos_Z)
		aa_kalman_compton_si_id.append(X_Y_flag)
		aa_kalman_compton_pos.append(Cluster_position)
		aa_kalman_compton_edep.append(energy_deposition)
		aa_kalman_compton_strip_number.append(number_of_strips_composing_the_cluster)
		aa_kalman_compton_pair.append(pair_flag)
		
		
		table_aa_kalman_cluster_compton.close()

		filenamedat_aa_kalman_ray = filepath+sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'
		
		table_aa_kalman_cluster_ray = open(filenamedat_aa_kalman_ray)
		
		lst = []
		
		for line in table_aa_kalman_cluster_ray:
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
		
		aa_kalman_ray_event_id.append(event_ID)
		aa_kalman_ray_theta_in.append(theta_input)
		aa_kalman_ray_phi_in.append(phi_input)
		aa_kalman_ray_ene_in.append(energy_input)
		aa_kalman_ray_plane_id.append(plane_ID)
		aa_kalman_ray_zpos.append(Pos_Z)
		aa_kalman_ray_si_id.append(X_Y_flag)
		aa_kalman_ray_pos.append(Cluster_position)
		aa_kalman_ray_edep.append(energy_deposition)
		aa_kalman_ray_strip_number.append(number_of_strips_composing_the_cluster)
		aa_kalman_ray_pair.append(pair_flag)
		
		
		table_aa_kalman_cluster_ray.close()

		filenamedat_S1_trk = filepath+sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'_TRACKER.dat'

		table_S1_trk = open(filenamedat_S1_trk)
		
		lst = []
		
		for line in table_S1_trk:
			lst += [line.split()]
		
		event_ID = [x[0] for x in lst]
		type = [x[1] for x in lst]
		Edep = [x[2] for x in lst]
		VolumeID = [x[3] for x in lst]
		X_pos = [x[4] for x in lst]
		Y_pos = [x[5] for x in lst]
		Z_pos = [x[6] for x in lst]
		
		S1_event_id_trk.append(event_ID)
		S1_type_trk.append(type)
		S1_e_dep_trk.append(Edep)
		S1_volume_id_trk.append(VolumeID)
		S1_x_trk.append(X_pos)
		S1_y_trk.append(Y_pos)
		S1_z_trk.append(Z_pos)
		
		
		table_S1_trk.close()

		filenamefits_raw = fits.open(filepath+'G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_raw[1].data

		rawData_event_id_temp = tbdata.field('EVT_ID')
		rawData_vol_id_temp = tbdata.field('VOL_ID')
		rawData_moth_id_temp = tbdata.field('MOTH_ID')
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
		rawData_trk_id_temp = tbdata.field('TRK_ID')
		rawData_part_id_temp = tbdata.field('PART_ID')


		rawData_event_id.append(rawData_event_id_temp)
		rawData_vol_id.append(rawData_vol_id_temp)
		rawData_moth_id.append(rawData_moth_id_temp)
		rawData_tray_id.append(rawData_tray_id_temp)
		rawData_plane_id.append(rawData_plane_id_temp)
		rawData_Strip_id_x.append(rawData_Strip_id_x_temp)
		rawData_Strip_id_y.append(rawData_Strip_id_y_temp)
		rawData_energy_dep.append(rawData_energy_dep_temp)
		rawData_ent_x.append(rawData_ent_x_temp)
		rawData_ent_y.append(rawData_ent_y_temp)
		rawData_ent_z.append(rawData_ent_z_temp)
		rawData_exit_x.append(rawData_exit_x_temp)
		rawData_exit_y.append(rawData_exit_y_temp)
		rawData_exit_z.append(rawData_exit_z_temp)
		rawData_child_id.append(rawData_child_id_temp)
		rawData_proc_id.append(rawData_proc_id_temp)
		rawData_trk_id.append(rawData_trk_id_temp)
		rawData_part_id.append(rawData_part_id_temp)


		filenamefits_raw.close()
		

		filenamefits_l0 = fits.open(filepath+'L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_l0[1].data

		L0TRACKER_Glob_event_id_temp = tbdata.field('EVT_ID')
		L0TRACKER_Glob_vol_id_temp = tbdata.field('VOL_ID')		
		L0TRACKER_Glob_moth_id_temp = tbdata.field('MOTH_ID')		
		L0TRACKER_Glob_tray_id_temp = tbdata.field('TRAY_ID')
		L0TRACKER_Glob_plane_id_temp = tbdata.field('PLANE_ID')
		L0TRACKER_Glob_Si_id_temp = tbdata.field('TRK_FLAG')
		L0TRACKER_Glob_Strip_id_temp = tbdata.field('STRIP_ID')
		L0TRACKER_Glob_pos_temp = tbdata.field('POS')
		L0TRACKER_Glob_zpos_temp = tbdata.field('ZPOS')
		L0TRACKER_Glob_energy_dep_temp = tbdata.field('E_DEP')
		L0TRACKER_Glob_pair_flag_temp = tbdata.field('PAIR_FLAG')
		

		L0TRACKER_Glob_event_id.append(L0TRACKER_Glob_event_id_temp)
		L0TRACKER_Glob_vol_id.append(L0TRACKER_Glob_vol_id_temp)
		L0TRACKER_Glob_moth_id.append(L0TRACKER_Glob_moth_id_temp)		
		L0TRACKER_Glob_tray_id.append(L0TRACKER_Glob_tray_id_temp)
		L0TRACKER_Glob_plane_id.append(L0TRACKER_Glob_plane_id_temp)
		L0TRACKER_Glob_Si_id.append(L0TRACKER_Glob_Si_id_temp)
		L0TRACKER_Glob_Strip_id.append(L0TRACKER_Glob_Strip_id_temp)
		L0TRACKER_Glob_pos.append(L0TRACKER_Glob_pos_temp)
		L0TRACKER_Glob_zpos.append(L0TRACKER_Glob_zpos_temp)
		L0TRACKER_Glob_energy_dep.append(L0TRACKER_Glob_energy_dep_temp)
		L0TRACKER_Glob_pair_flag.append(L0TRACKER_Glob_pair_flag_temp)

		filenamefits_l0.close()


		filenamefits_l05 = fits.open(filepath+'L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_l05[1].data

		L05TRACKER_Glob_event_id_cluster_temp = tbdata.field('EVT_ID')
		L05TRACKER_Glob_tray_id_cluster_temp = tbdata.field('TRAY_ID')
		L05TRACKER_Glob_plane_id_cluster_temp = tbdata.field('PLANE_ID')
		L05TRACKER_Glob_Si_id_cluster_temp = tbdata.field('TRK_FLAG')
		L05TRACKER_Glob_pos_cluster_temp = tbdata.field('POS')
		L05TRACKER_Glob_zpos_cluster_temp = tbdata.field('ZPOS')
		L05TRACKER_Glob_energy_dep_cluster_temp = tbdata.field('E_DEP')
		L05TRACKER_Glob_pair_flag_cluster_temp = tbdata.field('PAIR_FLAG')
		

		L05TRACKER_Glob_event_id_cluster.append(L05TRACKER_Glob_event_id_cluster_temp)
		L05TRACKER_Glob_tray_id_cluster.append(L05TRACKER_Glob_tray_id_cluster_temp)
		L05TRACKER_Glob_plane_id_cluster.append(L05TRACKER_Glob_plane_id_cluster_temp)
		L05TRACKER_Glob_Si_id_cluster.append(L05TRACKER_Glob_Si_id_cluster_temp)
		L05TRACKER_Glob_pos_cluster.append(L05TRACKER_Glob_pos_cluster_temp)
		L05TRACKER_Glob_zpos_cluster.append(L05TRACKER_Glob_zpos_cluster_temp)
		L05TRACKER_Glob_energy_cluster_dep.append(L05TRACKER_Glob_energy_dep_cluster_temp)
		L05TRACKER_Glob_pair_flag_cluster.append(L05TRACKER_Glob_pair_flag_cluster_temp)

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
		
		aa_fake_event_id.append(event_ID)
		aa_fake_theta_in.append(theta_input)
		aa_fake_phi_in.append(phi_input)
		aa_fake_ene_in.append(energy_input)
		aa_fake_plane_id.append(plane_ID)
		aa_fake_zpos.append(Pos_Z)
		aa_fake_si_id.append(X_Y_flag)
		aa_fake_pos.append(Cluster_position)
		aa_fake_edep.append(energy_deposition)
		aa_fake_strip_number.append(number_of_strips_composing_the_cluster)
		aa_fake_child_id.append(child_id)
		aa_fake_proc_id.append(proc_id)
		
		
		table_aa_fake.close()


	if cal_flag == 1:

		filenamefits_raw_cal = fits.open(filepath+'G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')

		tbdata = filenamefits_raw_cal[1].data

		rawData_event_id_cal_temp = tbdata.field('EVT_ID')
		rawData_vol_id_cal_temp = tbdata.field('VOL_ID')
		rawData_moth_id_cal_temp = tbdata.field('MOTH_ID')
		rawData_energy_dep_cal_temp = tbdata.field('E_DEP')
		rawData_ent_x_cal_temp = tbdata.field('X_ENT')
		rawData_ent_y_cal_temp = tbdata.field('Y_ENT')
		rawData_ent_z_cal_temp = tbdata.field('Z_ENT')
		rawData_exit_x_cal_temp = tbdata.field('X_EXIT')
		rawData_exit_y_cal_temp = tbdata.field('Y_EXIT')
		rawData_exit_z_cal_temp = tbdata.field('Z_EXIT')
		rawData_part_id_cal_temp = tbdata.field('PART_ID')
		rawData_trk_id_cal_temp = tbdata.field('TRK_ID')
		rawData_child_id_cal_temp = tbdata.field('CHILD_ID')
		rawData_proc_id_cal_temp = tbdata.field('PROC_ID')

		rawData_event_id_cal.append(rawData_event_id_cal_temp)
		rawData_vol_id_cal.append(rawData_vol_id_cal_temp)
		rawData_moth_id_cal.append(rawData_moth_id_cal_temp)
		rawData_energy_dep_cal.append(rawData_energy_dep_cal_temp)
		rawData_ent_x_cal.append(rawData_ent_x_cal_temp)
		rawData_ent_y_cal.append(rawData_ent_y_cal_temp)
		rawData_ent_z_cal.append(rawData_ent_z_cal_temp)
		rawData_exit_x_cal.append(rawData_exit_x_cal_temp)
		rawData_exit_y_cal.append(rawData_exit_y_cal_temp)
		rawData_exit_z_cal.append(rawData_exit_z_cal_temp)
		rawData_part_id_cal.append(rawData_part_id_cal_temp)
		rawData_trk_id_cal.append(rawData_trk_id_cal_temp)
		rawData_child_id_cal.append(rawData_child_id_cal_temp)
		rawData_proc_id_cal.append(rawData_proc_id_cal_temp)

		del tbdata
		filenamefits_raw_cal.close()
		filenamefits_cal = fits.open(filepath+'G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')

		tbdata = filenamefits_cal[1].data

		calInput_event_id_tot_temp = tbdata.field('EVT_ID')
		calInput_bar_id_tot_temp = tbdata.field('BAR_ID')
		calInput_bar_ene_tot_temp = tbdata.field('BAR_ENERGY')
		calInput_pair_flag_tot_temp = tbdata.field('PAIR_FLAG')


		calInput_event_id_tot.append(calInput_event_id_tot_temp)
		calInput_bar_id_tot.append(calInput_bar_id_tot_temp)
		calInput_bar_ene_tot.append(calInput_bar_ene_tot_temp)
		calInput_pair_flag_tot.append(calInput_pair_flag_tot_temp)

		del tbdata
		filenamefits_cal.close()
		
		if os.path.exists(filepath+'G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):


			filenamefits_cal_compton = fits.open(filepath+'G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
			tbdata = filenamefits_cal_compton[1].data

			calInput_event_id_tot_compton_temp = tbdata.field('EVT_ID')
			calInput_bar_id_tot_compton_temp = tbdata.field('BAR_ID')
			calInput_bar_ene_tot_compton_temp = tbdata.field('BAR_ENERGY')
			calInput_pair_flag_tot_compton_temp = tbdata.field('PAIR_FLAG')


			calInput_event_id_tot_compton.append(calInput_event_id_tot_compton_temp)
			calInput_bar_id_tot_compton.append(calInput_bar_id_tot_compton_temp)
			calInput_bar_ene_tot_compton.append(calInput_bar_ene_tot_compton_temp)
			calInput_pair_flag_tot_compton.append(calInput_pair_flag_tot_compton_temp)

			del tbdata
            filenamefits_cal_compton.close()

		else:
			pass




		if os.path.exists(filepath+'G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):


			filenamefits_cal_pair = fits.open(filepath+'G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
			tbdata = filenamefits_cal_pair[1].data

			calInput_event_id_tot_pair_temp = tbdata.field('EVT_ID')
			calInput_bar_id_tot_pair_temp = tbdata.field('BAR_ID')
			calInput_bar_ene_tot_pair_temp = tbdata.field('BAR_ENERGY')
			calInput_pair_flag_tot_pair_temp = tbdata.field('PAIR_FLAG')


			calInput_event_id_tot_pair.append(calInput_event_id_tot_pair_temp)
			calInput_bar_id_tot_pair.append(calInput_bar_id_tot_pair_temp)
			calInput_bar_ene_tot_pair.append(calInput_bar_ene_tot_pair_temp)
			calInput_pair_flag_tot_pair.append(calInput_pair_flag_tot_pair_temp)

			filenamefits_cal_pair.close()

		else:
			pass

		if os.path.exists(filepath+'G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):


			filenamefits_cal_ray = fits.open(filepath+'G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
			tbdata = filenamefits_cal_ray[1].data

			calInput_event_id_tot_ray_temp = tbdata.field('EVT_ID')
			calInput_bar_id_tot_ray_temp = tbdata.field('BAR_ID')
			calInput_bar_ene_tot_ray_temp = tbdata.field('BAR_ENERGY')
			calInput_pair_flag_tot_ray_temp = tbdata.field('PAIR_FLAG')


			calInput_event_id_tot_ray.append(calInput_event_id_tot_ray_temp)
			calInput_bar_id_tot_ray.append(calInput_bar_id_tot_ray_temp)
			calInput_bar_ene_tot_ray.append(calInput_bar_ene_tot_ray_temp)
			calInput_pair_flag_tot_ray.append(calInput_pair_flag_tot_ray_temp)

			filenamefits_cal_ray.close()

		else:
			pass

		filenamefits_cal_sum = fits.open(filepath+'SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_cal_sum[1].data

		calInputSum_event_id_tot_temp = tbdata.field('EVT_ID')
		calInputSum_bar_ene_tot_temp = tbdata.field('BAR_ENERGY')

		calInputSum_event_id_tot.append(calInputSum_event_id_tot_temp)
		calInputSum_bar_ene_tot.append(calInputSum_bar_ene_tot_temp)

		filenamefits_cal_sum.close()
		
		
		filenamedat_S1_cal = filepath+sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'_CAL.dat'

		table_S1_cal = open(filenamedat_S1_cal)
		
		lst = []
		
		for line in table_S1_cal:
			lst += [line.split()]
		
		event_ID = [x[0] for x in lst]
		type = [x[1] for x in lst]
		Edep = [x[2] for x in lst]
		VolumeID = [x[3] for x in lst]
		X_pos = [x[4] for x in lst]
		Y_pos = [x[5] for x in lst]
		Z_pos = [x[6] for x in lst]
		
		S1_event_id_cal.append(event_ID)
		S1_type_cal.append(type)
		S1_e_dep_cal.append(Edep)
		S1_volume_id_cal.append(VolumeID)
		S1_x_cal.append(X_pos)
		S1_y_cal.append(Y_pos)
		S1_z_cal.append(Z_pos)


	if ac_flag == 1:

		filenamefits_raw_ac = fits.open(filepath+'G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_raw_ac[1].data

		rawData_event_id_ac_temp = tbdata.field('EVT_ID')
		rawData_vol_id_ac_temp = tbdata.field('VOL_ID')
		rawData_moth_id_ac_temp = tbdata.field('MOTH_ID')
		rawData_energy_dep_ac_temp = tbdata.field('E_DEP')
		rawData_ent_x_ac_temp = tbdata.field('X_ENT')
		rawData_ent_y_ac_temp = tbdata.field('Y_ENT')
		rawData_ent_z_ac_temp = tbdata.field('Z_ENT')
		rawData_exit_x_ac_temp = tbdata.field('X_EXIT')
		rawData_exit_y_ac_temp = tbdata.field('Y_EXIT')
		rawData_exit_z_ac_temp = tbdata.field('Z_EXIT')
		rawData_part_id_ac_temp = tbdata.field('PART_ID')
		rawData_trk_id_ac_temp = tbdata.field('TRK_ID')
		rawData_child_id_ac_temp = tbdata.field('CHILD_ID')
		rawData_proc_id_ac_temp = tbdata.field('PROC_ID')
	
		rawData_event_id_ac.append(rawData_event_id_ac_temp)
		rawData_vol_id_ac.append(rawData_vol_id_ac_temp)
		rawData_moth_id_ac.append(rawData_moth_id_ac_temp)
		rawData_energy_dep_ac.append(rawData_energy_dep_ac_temp)
		rawData_ent_x_ac.append(rawData_ent_x_ac_temp)
		rawData_ent_y_ac.append(rawData_ent_y_ac_temp)
		rawData_ent_z_ac.append(rawData_ent_z_ac_temp)
		rawData_exit_x_ac.append(rawData_exit_x_ac_temp)
		rawData_exit_y_ac.append(rawData_exit_y_ac_temp)
		rawData_exit_z_ac.append(rawData_exit_z_ac_temp)
		rawData_part_id_ac.append(rawData_part_id_ac_temp)
		rawData_trk_id_ac.append(rawData_trk_id_ac_temp)
		rawData_child_id_ac.append(rawData_child_id_ac_temp)
		rawData_proc_id_ac.append(rawData_proc_id_ac_temp)

		filenamefits_raw_ac.close()

		
		filenamefits_ac = fits.open(filepath+'G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
		tbdata = filenamefits_ac[1].data

		acInput_event_id_tot_temp = tbdata.field('EVT_ID')
		acInput_AC_panel_temp = tbdata.field('AC_PANEL')
		acInput_AC_subpanel_temp = tbdata.field('AC_SUBPANEL')
		acInput_energy_dep_tot_temp = tbdata.field('E_DEP')
		acInput_pair_flag_tot_temp = tbdata.field('PAIR_FLAG')


		acInput_event_id_tot.append(acInput_event_id_tot_temp)
		acInput_AC_panel.append(acInput_AC_panel_temp)
		acInput_AC_subpanel.append(acInput_AC_subpanel_temp)
		acInput_energy_dep_tot.append(acInput_energy_dep_tot_temp)
		acInput_pair_flag_tot.append(acInput_pair_flag_tot_temp)

		filenamefits_ac.close()
		

		if os.path.exists(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):

			filenamefits_ac_compton = fits.open(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
			tbdata = filenamefits_ac_compton[1].data

			acInput_event_id_tot_compton_temp = tbdata.field('EVT_ID')
			acInput_AC_panel_compton_temp = tbdata.field('AC_PANEL')
			acInput_AC_subpanel_compton_temp = tbdata.field('AC_SUBPANEL')
			acInput_energy_dep_tot_compton_temp = tbdata.field('E_DEP')
			acInput_pair_flag_tot_compton_temp = tbdata.field('PAIR_FLAG')


			acInput_event_id_tot_compton.append(acInput_event_id_tot_compton_temp)
			acInput_AC_panel_compton.append(acInput_AC_panel_compton_temp)
			acInput_AC_subpanel_compton.append(acInput_AC_subpanel_compton_temp)
			acInput_energy_dep_tot_compton.append(acInput_energy_dep_tot_compton_temp)
			acInput_pair_flag_tot_compton.append(acInput_pair_flag_tot_compton_temp)

			filenamefits_ac_compton.close()		
		
		else:
			pass


		if os.path.exists(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):

			filenamefits_ac_pair = fits.open(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
			tbdata = filenamefits_ac_pair[1].data

			acInput_event_id_tot_pair_temp = tbdata.field('EVT_ID')
			acInput_AC_panel_pair_temp = tbdata.field('AC_PANEL')
			acInput_AC_subpanel_pair_temp = tbdata.field('AC_SUBPANEL')
			acInput_energy_dep_tot_pair_temp = tbdata.field('E_DEP')
			acInput_pair_flag_tot_pair_temp = tbdata.field('PAIR_FLAG')


			acInput_event_id_tot_pair.append(acInput_event_id_tot_pair_temp)
			acInput_AC_panel_pair.append(acInput_AC_panel_pair_temp)
			acInput_AC_subpanel_pair.append(acInput_AC_subpanel_pair_temp)
			acInput_energy_dep_tot_pair.append(acInput_energy_dep_tot_pair_temp)
			acInput_pair_flag_tot_pair.append(acInput_pair_flag_tot_pair_temp)

			filenamefits_ac_pair.close()		

		else:
			pass
			
		if os.path.exists(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):

			filenamefits_ac_ray = fits.open(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
   
			tbdata = filenamefits_ac_ray[1].data

			acInput_event_id_tot_ray_temp = tbdata.field('EVT_ID')
			acInput_AC_panel_ray_temp = tbdata.field('AC_PANEL')
			acInput_AC_subpanel_ray_temp = tbdata.field('AC_SUBPANEL')
			acInput_energy_dep_tot_ray_temp = tbdata.field('E_DEP')
			acInput_pair_flag_tot_ray_temp = tbdata.field('PAIR_FLAG')


			acInput_event_id_tot_ray.append(acInput_event_id_tot_ray_temp)
			acInput_AC_panel_ray.append(acInput_AC_panel_ray_temp)
			acInput_AC_subpanel_ray.append(acInput_AC_subpanel_ray_temp)
			acInput_energy_dep_tot_ray.append(acInput_energy_dep_tot_ray_temp)
			acInput_pair_flag_tot_ray.append(acInput_pair_flag_tot_ray_temp)

			filenamefits_ac_ray.close()		
		
		else:
			pass

		filenamedat_S1_ac = filepath+sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'_AC.dat'

		table_S1_ac = open(filenamedat_S1_ac)
		
		lst = []
		
		for line in table_S1_ac:
			lst += [line.split()]
		
		event_ID = [x[0] for x in lst]
		type = [x[1] for x in lst]
		Edep = [x[2] for x in lst]
		VolumeID = [x[3] for x in lst]
		X_pos = [x[4] for x in lst]
		Y_pos = [x[5] for x in lst]
		Z_pos = [x[6] for x in lst]
		
		S1_event_id_ac.append(event_ID)
		S1_type_ac.append(type)
		S1_e_dep_ac.append(Edep)
		S1_volume_id_ac.append(VolumeID)
		S1_x_ac.append(X_pos)
		S1_y_ac.append(Y_pos)
		S1_z_ac.append(Z_pos)

				
#############################################

	ifile = ifile + 1



aa_strip_event_id = np.ma.concatenate(aa_strip_event_id)
aa_strip_theta_in = np.ma.concatenate(aa_strip_theta_in)
aa_strip_phi_in = np.ma.concatenate(aa_strip_phi_in)
aa_strip_ene_in = np.ma.concatenate(aa_strip_ene_in)
aa_strip_plane_id = np.ma.concatenate(aa_strip_plane_id)
aa_strip_zpos = np.ma.concatenate(aa_strip_zpos)
aa_strip_si_id = np.ma.concatenate(aa_strip_si_id)
aa_strip_strip_id = np.ma.concatenate(aa_strip_strip_id)
aa_strip_pos = np.ma.concatenate(aa_strip_pos)
aa_strip_edep = np.ma.concatenate(aa_strip_edep)
aa_strip_pair = np.ma.concatenate(aa_strip_pair)

aa_strip_event_id = np.array((aa_strip_event_id), dtype=np.int64)
aa_strip_theta_in = np.array((aa_strip_theta_in), dtype=np.int64)
aa_strip_phi_in = np.array((aa_strip_phi_in), dtype=np.int64)
aa_strip_ene_in = np.array(aa_strip_ene_in)
aa_strip_plane_id = np.array((aa_strip_plane_id), dtype=np.int64)
aa_strip_zpos = np.array((aa_strip_zpos), dtype=np.float64)
aa_strip_si_id = np.array((aa_strip_si_id), dtype=np.int64)
aa_strip_strip_id = np.array((aa_strip_strip_id), dtype=np.int64)
aa_strip_pos = np.array((aa_strip_pos), dtype=np.float64)
aa_strip_edep = np.array((aa_strip_edep), dtype=np.float64)
aa_strip_pair = np.array((aa_strip_pair), dtype=np.int64)


aa_kalman_event_id = np.ma.concatenate(aa_kalman_event_id)
aa_kalman_theta_in = np.ma.concatenate(aa_kalman_theta_in)
aa_kalman_phi_in = np.ma.concatenate(aa_kalman_phi_in)
aa_kalman_ene_in = np.ma.concatenate(aa_kalman_ene_in)
aa_kalman_plane_id = np.ma.concatenate(aa_kalman_plane_id)
aa_kalman_zpos = np.ma.concatenate(aa_kalman_zpos)
aa_kalman_si_id = np.ma.concatenate(aa_kalman_si_id)
aa_kalman_pos = np.ma.concatenate(aa_kalman_pos)
aa_kalman_edep = np.ma.concatenate(aa_kalman_edep)
aa_kalman_strip_number = np.ma.concatenate(aa_kalman_strip_number)
aa_kalman_pair = np.ma.concatenate(aa_kalman_pair)

aa_kalman_event_id = np.array((aa_kalman_event_id), dtype=np.int64)
aa_kalman_theta_in = np.array((aa_kalman_theta_in), dtype=np.int64)
aa_kalman_phi_in = np.array((aa_kalman_phi_in), dtype=np.int64)
aa_kalman_ene_in = np.array(aa_kalman_ene_in)
aa_kalman_plane_id = np.array((aa_kalman_plane_id), dtype=np.int64)
aa_kalman_zpos = np.array((aa_kalman_zpos), dtype=np.float64)
aa_kalman_si_id = np.array((aa_kalman_si_id), dtype=np.int64)
aa_kalman_pos = np.array((aa_kalman_pos), dtype=np.float64)
aa_kalman_edep = np.array((aa_kalman_edep), dtype=np.float64)
aa_kalman_strip_number = np.array((aa_kalman_strip_number), dtype=np.int64)
aa_kalman_pair = np.array((aa_kalman_pair), dtype=np.int64)


aa_kalman_pair_event_id = np.ma.concatenate(aa_kalman_pair_event_id)
aa_kalman_pair_theta_in = np.ma.concatenate(aa_kalman_pair_theta_in)
aa_kalman_pair_phi_in = np.ma.concatenate(aa_kalman_pair_phi_in)
aa_kalman_pair_ene_in = np.ma.concatenate(aa_kalman_pair_ene_in)
aa_kalman_pair_plane_id = np.ma.concatenate(aa_kalman_pair_plane_id)
aa_kalman_pair_zpos = np.ma.concatenate(aa_kalman_pair_zpos)
aa_kalman_pair_si_id = np.ma.concatenate(aa_kalman_pair_si_id)
aa_kalman_pair_pos = np.ma.concatenate(aa_kalman_pair_pos)
aa_kalman_pair_edep = np.ma.concatenate(aa_kalman_pair_edep)
aa_kalman_pair_strip_number = np.ma.concatenate(aa_kalman_pair_strip_number)
aa_kalman_pair_pair = np.ma.concatenate(aa_kalman_pair_pair)

aa_kalman_pair_event_id = np.array((aa_kalman_pair_event_id), dtype=np.int64)
aa_kalman_pair_theta_in = np.array((aa_kalman_pair_theta_in), dtype=np.int64)
aa_kalman_pair_phi_in = np.array((aa_kalman_pair_phi_in), dtype=np.int64)
aa_kalman_pair_ene_in = np.array(aa_kalman_pair_ene_in)
aa_kalman_pair_plane_id = np.array((aa_kalman_pair_plane_id), dtype=np.int64)
aa_kalman_pair_zpos = np.array((aa_kalman_pair_zpos), dtype=np.float64)
aa_kalman_pair_si_id = np.array((aa_kalman_pair_si_id), dtype=np.int64)
aa_kalman_pair_pos = np.array((aa_kalman_pair_pos), dtype=np.float64)
aa_kalman_pair_edep = np.array((aa_kalman_pair_edep), dtype=np.float64)
aa_kalman_pair_strip_number = np.array((aa_kalman_pair_strip_number), dtype=np.int64)
aa_kalman_pair_pair = np.array((aa_kalman_pair_pair), dtype=np.int64)


aa_kalman_compton_event_id = np.ma.concatenate(aa_kalman_compton_event_id)
aa_kalman_compton_theta_in = np.ma.concatenate(aa_kalman_compton_theta_in)
aa_kalman_compton_phi_in = np.ma.concatenate(aa_kalman_compton_phi_in)
aa_kalman_compton_ene_in = np.ma.concatenate(aa_kalman_compton_ene_in)
aa_kalman_compton_plane_id = np.ma.concatenate(aa_kalman_compton_plane_id)
aa_kalman_compton_zpos = np.ma.concatenate(aa_kalman_compton_zpos)
aa_kalman_compton_si_id = np.ma.concatenate(aa_kalman_compton_si_id)
aa_kalman_compton_pos = np.ma.concatenate(aa_kalman_compton_pos)
aa_kalman_compton_edep = np.ma.concatenate(aa_kalman_compton_edep)
aa_kalman_compton_strip_number = np.ma.concatenate(aa_kalman_compton_strip_number)
aa_kalman_compton_pair = np.ma.concatenate(aa_kalman_compton_pair)

aa_kalman_compton_event_id = np.array((aa_kalman_compton_event_id), dtype=np.int64)
aa_kalman_compton_theta_in = np.array((aa_kalman_compton_theta_in), dtype=np.int64)
aa_kalman_compton_phi_in = np.array((aa_kalman_compton_phi_in), dtype=np.int64)
aa_kalman_compton_ene_in = np.array(aa_kalman_compton_ene_in)
aa_kalman_compton_plane_id = np.array((aa_kalman_compton_plane_id), dtype=np.int64)
aa_kalman_compton_zpos = np.array((aa_kalman_compton_zpos), dtype=np.float64)
aa_kalman_compton_si_id = np.array((aa_kalman_compton_si_id), dtype=np.int64)
aa_kalman_compton_pos = np.array((aa_kalman_compton_pos), dtype=np.float64)
aa_kalman_compton_edep = np.array((aa_kalman_compton_edep), dtype=np.float64)
aa_kalman_compton_strip_number = np.array((aa_kalman_compton_strip_number), dtype=np.int64)
aa_kalman_compton_pair = np.array((aa_kalman_compton_pair), dtype=np.int64)

aa_kalman_ray_event_id = np.ma.concatenate(aa_kalman_ray_event_id)
aa_kalman_ray_theta_in = np.ma.concatenate(aa_kalman_ray_theta_in)
aa_kalman_ray_phi_in = np.ma.concatenate(aa_kalman_ray_phi_in)
aa_kalman_ray_ene_in = np.ma.concatenate(aa_kalman_ray_ene_in)
aa_kalman_ray_plane_id = np.ma.concatenate(aa_kalman_ray_plane_id)
aa_kalman_ray_zpos = np.ma.concatenate(aa_kalman_ray_zpos)
aa_kalman_ray_si_id = np.ma.concatenate(aa_kalman_ray_si_id)
aa_kalman_ray_pos = np.ma.concatenate(aa_kalman_ray_pos)
aa_kalman_ray_edep = np.ma.concatenate(aa_kalman_ray_edep)
aa_kalman_ray_strip_number = np.ma.concatenate(aa_kalman_ray_strip_number)
aa_kalman_ray_pair = np.ma.concatenate(aa_kalman_ray_pair)

aa_kalman_ray_event_id = np.array((aa_kalman_ray_event_id), dtype=np.int64)
aa_kalman_ray_theta_in = np.array((aa_kalman_ray_theta_in), dtype=np.int64)
aa_kalman_ray_phi_in = np.array((aa_kalman_ray_phi_in), dtype=np.int64)
aa_kalman_ray_ene_in = np.array(aa_kalman_ray_ene_in)
aa_kalman_ray_plane_id = np.array((aa_kalman_ray_plane_id), dtype=np.int64)
aa_kalman_ray_zpos = np.array((aa_kalman_ray_zpos), dtype=np.float64)
aa_kalman_ray_si_id = np.array((aa_kalman_ray_si_id), dtype=np.int64)
aa_kalman_ray_pos = np.array((aa_kalman_ray_pos), dtype=np.float64)
aa_kalman_ray_edep = np.array((aa_kalman_ray_edep), dtype=np.float64)
aa_kalman_ray_strip_number = np.array((aa_kalman_ray_strip_number), dtype=np.int64)
aa_kalman_ray_pair = np.array((aa_kalman_ray_pair), dtype=np.int64)

S1_event_id_trk = np.ma.concatenate(S1_event_id_trk)
S1_type_trk = np.ma.concatenate(S1_type_trk)
S1_e_dep_trk = np.ma.concatenate(S1_e_dep_trk)
S1_volume_id_trk = np.ma.concatenate(S1_volume_id_trk)
S1_x_trk = np.ma.concatenate(S1_x_trk)
S1_y_trk = np.ma.concatenate(S1_y_trk)
S1_z_trk = np.ma.concatenate(S1_z_trk)

S1_event_id_trk = np.array((S1_event_id_trk), dtype=np.int64)
S1_type_trk = np.array((S1_type_trk), dtype=np.int64)
S1_e_dep_trk = np.array((S1_e_dep_trk), dtype=np.float64)
S1_volume_id_trk = np.array((S1_volume_id_trk), dtype=np.int64)
S1_x_trk = np.array((S1_x_trk), dtype=np.float64)
S1_y_trk = np.array((S1_y_trk), dtype=np.float64)
S1_z_trk = np.array((S1_z_trk), dtype=np.float64)

rawData_event_id = np.ma.concatenate(rawData_event_id)
rawData_tray_id = np.ma.concatenate(rawData_tray_id)
rawData_plane_id = np.ma.concatenate(rawData_plane_id)
rawData_Strip_id_x = np.ma.concatenate(rawData_Strip_id_x)
rawData_Strip_id_y = np.ma.concatenate(rawData_Strip_id_y)
rawData_energy_dep = np.ma.concatenate(rawData_energy_dep)
rawData_ent_x = np.ma.concatenate(rawData_ent_x)
rawData_ent_y = np.ma.concatenate(rawData_ent_y)
rawData_ent_z = np.ma.concatenate(rawData_ent_z)
rawData_exit_x = np.ma.concatenate(rawData_exit_x)
rawData_exit_y = np.ma.concatenate(rawData_exit_y)
rawData_exit_z = np.ma.concatenate(rawData_exit_z)
rawData_child_id = np.ma.concatenate(rawData_child_id)
rawData_proc_id = np.ma.concatenate(rawData_proc_id)
rawData_trk_id = np.ma.concatenate(rawData_trk_id)
rawData_vol_id = np.ma.concatenate(rawData_vol_id)
rawData_moth_id = np.ma.concatenate(rawData_moth_id)
rawData_part_id = np.ma.concatenate(rawData_part_id)

rawData_event_id = np.array(rawData_event_id)
rawData_tray_id = np.array(rawData_tray_id)
rawData_plane_id = np.array(rawData_plane_id)
rawData_Strip_id_x = np.array(rawData_Strip_id_x)
rawData_Strip_id_y = np.array(rawData_Strip_id_y)
rawData_energy_dep = np.array(rawData_energy_dep)
rawData_ent_x = np.array(rawData_ent_x)
rawData_ent_y = np.array(rawData_ent_y)
rawData_ent_z = np.array(rawData_ent_z)
rawData_exit_x = np.array(rawData_exit_x)
rawData_exit_y = np.array(rawData_exit_y)
rawData_exit_z = np.array(rawData_exit_z)
rawData_child_id = np.array(rawData_child_id)
rawData_proc_id = np.array(rawData_proc_id)
rawData_trk_id = np.array(rawData_trk_id)
rawData_vol_id = np.array(rawData_vol_id)
rawData_moth_id = np.array(rawData_moth_id)
rawData_part_id = np.array(rawData_part_id)


L0TRACKER_Glob_event_id = np.ma.concatenate(L0TRACKER_Glob_event_id)
L0TRACKER_Glob_vol_id = np.ma.concatenate(L0TRACKER_Glob_vol_id)		
L0TRACKER_Glob_moth_id = np.ma.concatenate(L0TRACKER_Glob_moth_id)		
L0TRACKER_Glob_tray_id = np.ma.concatenate(L0TRACKER_Glob_tray_id)
L0TRACKER_Glob_plane_id = np.ma.concatenate(L0TRACKER_Glob_plane_id)
L0TRACKER_Glob_Si_id = np.ma.concatenate(L0TRACKER_Glob_Si_id)
L0TRACKER_Glob_Strip_id = np.ma.concatenate(L0TRACKER_Glob_Strip_id)
L0TRACKER_Glob_pos = np.ma.concatenate(L0TRACKER_Glob_pos)
L0TRACKER_Glob_zpos = np.ma.concatenate(L0TRACKER_Glob_zpos)
L0TRACKER_Glob_energy_dep = np.ma.concatenate(L0TRACKER_Glob_energy_dep)
L0TRACKER_Glob_pair_flag = np.ma.concatenate(L0TRACKER_Glob_pair_flag)

L0TRACKER_Glob_event_id = np.array(L0TRACKER_Glob_event_id)
L0TRACKER_Glob_vol_id = np.array(L0TRACKER_Glob_vol_id)		
L0TRACKER_Glob_moth_id = np.array(L0TRACKER_Glob_moth_id)		
L0TRACKER_Glob_tray_id = np.array(L0TRACKER_Glob_tray_id)
L0TRACKER_Glob_plane_id = np.array(L0TRACKER_Glob_plane_id)
L0TRACKER_Glob_Si_id = np.array(L0TRACKER_Glob_Si_id)
L0TRACKER_Glob_Strip_id = np.array(L0TRACKER_Glob_Strip_id)
L0TRACKER_Glob_pos = np.array(L0TRACKER_Glob_pos)
L0TRACKER_Glob_zpos = np.array(L0TRACKER_Glob_zpos)
L0TRACKER_Glob_energy_dep = np.array(L0TRACKER_Glob_energy_dep)
L0TRACKER_Glob_pair_flag = np.array(L0TRACKER_Glob_pair_flag)


L05TRACKER_Glob_event_id_cluster = np.ma.concatenate(L05TRACKER_Glob_event_id_cluster)
L05TRACKER_Glob_tray_id_cluster = np.ma.concatenate(L05TRACKER_Glob_tray_id_cluster)
L05TRACKER_Glob_plane_id_cluster = np.ma.concatenate(L05TRACKER_Glob_plane_id_cluster)
L05TRACKER_Glob_Si_id_cluster = np.ma.concatenate(L05TRACKER_Glob_Si_id_cluster)
L05TRACKER_Glob_pos_cluster = np.ma.concatenate(L05TRACKER_Glob_pos_cluster)
L05TRACKER_Glob_zpos_cluster = np.ma.concatenate(L05TRACKER_Glob_zpos_cluster)
L05TRACKER_Glob_energy_cluster_dep = np.ma.concatenate(L05TRACKER_Glob_energy_cluster_dep)
L05TRACKER_Glob_pair_flag_cluster = np.ma.concatenate(L05TRACKER_Glob_pair_flag_cluster)

L05TRACKER_Glob_event_id_cluster = np.array(L05TRACKER_Glob_event_id_cluster)
L05TRACKER_Glob_tray_id_cluster = np.array(L05TRACKER_Glob_tray_id_cluster)
L05TRACKER_Glob_plane_id_cluster = np.array(L05TRACKER_Glob_plane_id_cluster)
L05TRACKER_Glob_Si_id_cluster = np.array(L05TRACKER_Glob_Si_id_cluster)
L05TRACKER_Glob_pos_cluster = np.array(L05TRACKER_Glob_pos_cluster)
L05TRACKER_Glob_zpos_cluster = np.array(L05TRACKER_Glob_zpos_cluster)
L05TRACKER_Glob_energy_cluster_dep = np.array(L05TRACKER_Glob_energy_cluster_dep)
L05TRACKER_Glob_pair_flag_cluster = np.array(L05TRACKER_Glob_pair_flag_cluster)

if isStrip == 0:

	aa_fake_event_id = np.ma.concatenate(aa_fake_event_id)
	aa_fake_theta_in = np.ma.concatenate(aa_fake_theta_in)
	aa_fake_phi_in = np.ma.concatenate(aa_fake_phi_in)
	aa_fake_ene_in = np.ma.concatenate(aa_fake_ene_in)
	aa_fake_plane_id = np.ma.concatenate(aa_fake_plane_id)
	aa_fake_zpos = np.ma.concatenate(aa_fake_zpos)
	aa_fake_si_id = np.ma.concatenate(aa_fake_si_id)
	aa_fake_pos = np.ma.concatenate(aa_fake_pos)
	aa_fake_edep = np.ma.concatenate(aa_fake_edep)
	aa_fake_strip_number = np.ma.concatenate(aa_fake_strip_number)
	aa_fake_child_id = np.ma.concatenate(aa_fake_child_id)
	aa_fake_proc_id = np.ma.concatenate(aa_fake_proc_id)

	aa_fake_event_id = np.array((aa_fake_event_id), dtype=np.int64)
	aa_fake_theta_in = np.array(aa_fake_theta_in)
	aa_fake_phi_in = np.array(aa_fake_phi_in)
	aa_fake_ene_in = np.array(aa_fake_ene_in)
	aa_fake_plane_id = np.array(aa_fake_plane_id)
	aa_fake_zpos = np.array(aa_fake_zpos)
	aa_fake_si_id = np.array(aa_fake_si_id)
	aa_fake_pos = np.array(aa_fake_pos)
	aa_fake_edep = np.array(aa_fake_edep)
	aa_fake_strip_number = np.array(aa_fake_strip_number)
	aa_fake_child_id = np.array(aa_fake_child_id)
	aa_fake_proc_id = np.array(aa_fake_proc_id)


if cal_flag == 1:

	rawData_event_id_cal = np.ma.concatenate(rawData_event_id_cal)
	rawData_energy_dep_cal = np.ma.concatenate(rawData_energy_dep_cal)
	rawData_ent_x_cal = np.ma.concatenate(rawData_ent_x_cal)
	rawData_ent_y_cal = np.ma.concatenate(rawData_ent_y_cal)
	rawData_ent_z_cal = np.ma.concatenate(rawData_ent_z_cal)
	rawData_exit_x_cal = np.ma.concatenate(rawData_exit_x_cal)
	rawData_exit_y_cal = np.ma.concatenate(rawData_exit_y_cal)
	rawData_exit_z_cal = np.ma.concatenate(rawData_exit_z_cal)
	rawData_child_id_cal = np.ma.concatenate(rawData_child_id_cal)
	rawData_proc_id_cal = np.ma.concatenate(rawData_proc_id_cal)
	rawData_trk_id_cal = np.ma.concatenate(rawData_trk_id_cal)
	rawData_vol_id_cal = np.ma.concatenate(rawData_vol_id_cal)
	rawData_moth_id_cal = np.ma.concatenate(rawData_moth_id_cal)
	rawData_part_id_cal = np.ma.concatenate(rawData_part_id_cal)
	
	rawData_event_id_cal = np.array(rawData_event_id_cal)
	rawData_energy_dep_cal = np.array(rawData_energy_dep_cal)
	rawData_ent_x_cal = np.array(rawData_ent_x_cal)
	rawData_ent_y_cal = np.array(rawData_ent_y_cal)
	rawData_ent_z_cal = np.array(rawData_ent_z_cal)
	rawData_exit_x_cal = np.array(rawData_exit_x_cal)
	rawData_exit_y_cal = np.array(rawData_exit_y_cal)
	rawData_exit_z_cal = np.array(rawData_exit_z_cal)
	rawData_child_id_cal = np.array(rawData_child_id_cal)
	rawData_proc_id_cal = np.array(rawData_proc_id_cal)
	rawData_trk_id_cal = np.array(rawData_trk_id_cal)
	rawData_vol_id_cal = np.array(rawData_vol_id_cal)
	rawData_moth_id_cal = np.array(rawData_moth_id_cal)
	rawData_part_id_cal = np.array(rawData_part_id_cal)

	
	calInput_event_id_tot = np.ma.concatenate(calInput_event_id_tot)
	calInput_bar_id_tot = np.ma.concatenate(calInput_bar_id_tot)
	calInput_bar_ene_tot = np.ma.concatenate(calInput_bar_ene_tot)
	calInput_pair_flag_tot = np.ma.concatenate(calInput_pair_flag_tot)	
	
	calInput_event_id_tot = np.array(calInput_event_id_tot)
	calInput_bar_id_tot = np.array(calInput_bar_id_tot)
	calInput_bar_ene_tot = np.array(calInput_bar_ene_tot)
	calInput_pair_flag_tot = np.array(calInput_pair_flag_tot)
	
	if (len(calInput_event_id_tot_compton) > 0):
		calInput_event_id_tot_compton = np.ma.concatenate(calInput_event_id_tot_compton)
		calInput_bar_id_tot_compton = np.ma.concatenate(calInput_bar_id_tot_compton)
		calInput_bar_ene_tot_compton = np.ma.concatenate(calInput_bar_ene_tot_compton)
		calInput_pair_flag_tot_compton = np.ma.concatenate(calInput_pair_flag_tot_compton)	
	
	calInput_event_id_tot_compton = np.array(calInput_event_id_tot_compton)
	calInput_bar_id_tot_compton = np.array(calInput_bar_id_tot_compton)
	calInput_bar_ene_tot_compton = np.array(calInput_bar_ene_tot_compton)
	calInput_pair_flag_tot_compton = np.array(calInput_pair_flag_tot_compton)

	if (len(calInput_event_id_tot_pair) > 0):
		calInput_event_id_tot_pair = np.ma.concatenate(calInput_event_id_tot_pair)
		calInput_bar_id_tot_pair = np.ma.concatenate(calInput_bar_id_tot_pair)
		calInput_bar_ene_tot_pair = np.ma.concatenate(calInput_bar_ene_tot_pair)
		calInput_pair_flag_tot_pair = np.ma.concatenate(calInput_pair_flag_tot_pair)	
	
	calInput_event_id_tot_pair = np.array(calInput_event_id_tot_pair)
	calInput_bar_id_tot_pair = np.array(calInput_bar_id_tot_pair)
	calInput_bar_ene_tot_pair = np.array(calInput_bar_ene_tot_pair)
	calInput_pair_flag_tot_pair = np.array(calInput_pair_flag_tot_pair)

	if (len(calInput_event_id_tot_ray) > 0):
		calInput_event_id_tot_ray = np.ma.concatenate(calInput_event_id_tot_ray)
		calInput_bar_id_tot_ray = np.ma.concatenate(calInput_bar_id_tot_ray)
		calInput_bar_ene_tot_ray = np.ma.concatenate(calInput_bar_ene_tot_ray)
		calInput_pair_flag_tot_ray = np.ma.concatenate(calInput_pair_flag_tot_ray)	
	
	calInput_event_id_tot_ray = np.array(calInput_event_id_tot_ray)
	calInput_bar_id_tot_ray = np.array(calInput_bar_id_tot_ray)
	calInput_bar_ene_tot_ray = np.array(calInput_bar_ene_tot_ray)
	calInput_pair_flag_tot_ray = np.array(calInput_pair_flag_tot_ray)

	
	calInputSum_event_id_tot = np.ma.concatenate(calInputSum_event_id_tot)
	calInputSum_bar_ene_tot = np.ma.concatenate(calInputSum_bar_ene_tot)
	
	calInputSum_event_id_tot = np.array(calInputSum_event_id_tot)
	calInputSum_bar_ene_tot = np.array(calInputSum_bar_ene_tot)
	
	
	
	S1_event_id_cal = np.ma.concatenate(S1_event_id_cal)
	S1_type_cal = np.ma.concatenate(S1_type_cal)
	S1_e_dep_cal = np.ma.concatenate(S1_e_dep_cal)
	S1_volume_id_cal = np.ma.concatenate(S1_volume_id_cal)
	S1_x_cal = np.ma.concatenate(S1_x_cal)
	S1_y_cal = np.ma.concatenate(S1_y_cal)
	S1_z_cal = np.ma.concatenate(S1_z_cal)
	
	S1_event_id_cal = np.array((S1_event_id_cal), dtype=np.int64)
	S1_type_cal = np.array((S1_type_cal), dtype=np.int64)
	S1_e_dep_cal = np.array((S1_e_dep_cal), dtype=np.float64)
	S1_volume_id_cal = np.array((S1_volume_id_cal), dtype=np.int64)
	S1_x_cal = np.array((S1_x_cal), dtype=np.float64)
	S1_y_cal = np.array((S1_y_cal), dtype=np.float64)
	S1_z_cal = np.array((S1_z_cal), dtype=np.float64)

if ac_flag == 1:

	
	rawData_event_id_ac = np.ma.concatenate(rawData_event_id_ac)
	rawData_energy_dep_ac = np.ma.concatenate(rawData_energy_dep_ac)
	rawData_ent_x_ac = np.ma.concatenate(rawData_ent_x_ac)
	rawData_ent_y_ac = np.ma.concatenate(rawData_ent_y_ac)
	rawData_ent_z_ac = np.ma.concatenate(rawData_ent_z_ac)
	rawData_exit_x_ac = np.ma.concatenate(rawData_exit_x_ac)
	rawData_exit_y_ac = np.ma.concatenate(rawData_exit_y_ac)
	rawData_exit_z_ac = np.ma.concatenate(rawData_exit_z_ac)
	rawData_child_id_ac = np.ma.concatenate(rawData_child_id_ac)
	rawData_proc_id_ac = np.ma.concatenate(rawData_proc_id_ac)
	rawData_trk_id_ac = np.ma.concatenate(rawData_trk_id_ac)
	rawData_vol_id_ac = np.ma.concatenate(rawData_vol_id_ac)
	rawData_moth_id_ac = np.ma.concatenate(rawData_moth_id_ac)
	rawData_part_id_ac = np.ma.concatenate(rawData_part_id_ac)
	

	rawData_event_id_ac = np.array(rawData_event_id_ac)
	rawData_energy_dep_ac = np.array(rawData_energy_dep_ac)
	rawData_ent_x_ac = np.array(rawData_ent_x_ac)
	rawData_ent_y_ac = np.array(rawData_ent_y_ac)
	rawData_ent_z_ac = np.array(rawData_ent_z_ac)
	rawData_exit_x_ac = np.array(rawData_exit_x_ac)
	rawData_exit_y_ac = np.array(rawData_exit_y_ac)
	rawData_exit_z_ac = np.array(rawData_exit_z_ac)
	rawData_child_id_ac = np.array(rawData_child_id_ac)
	rawData_proc_id_ac = np.array(rawData_proc_id_ac)
	rawData_trk_id_ac = np.array(rawData_trk_id_ac)
	rawData_vol_id_ac = np.array(rawData_vol_id_ac)
	rawData_moth_id_ac = np.array(rawData_moth_id_ac)
	rawData_part_id_ac = np.array(rawData_part_id_ac)



	
	acInput_event_id_tot = np.ma.concatenate(acInput_event_id_tot)
	acInput_AC_panel = np.ma.concatenate(acInput_AC_panel)
	acInput_AC_subpanel = np.ma.concatenate(acInput_AC_subpanel)
	acInput_energy_dep_tot = np.ma.concatenate(acInput_energy_dep_tot)
	acInput_pair_flag_tot = np.ma.concatenate(acInput_pair_flag_tot)
	
		
	acInput_event_id_tot = np.array(acInput_event_id_tot)
	acInput_AC_panel = np.array(acInput_AC_panel)
	acInput_AC_subpanel = np.array(acInput_AC_subpanel)
	acInput_energy_dep_tot = np.array(acInput_energy_dep_tot)
	acInput_pair_flag_tot = np.array(acInput_pair_flag_tot)


	if (len(acInput_event_id_tot_compton) > 0):
		acInput_event_id_tot_compton = np.ma.concatenate(acInput_event_id_tot_compton)
		acInput_AC_panel_compton = np.ma.concatenate(acInput_AC_panel_compton)
		acInput_AC_subpanel_compton = np.ma.concatenate(acInput_AC_subpanel_compton)
		acInput_energy_dep_tot_compton = np.ma.concatenate(acInput_energy_dep_tot_compton)
		acInput_pair_flag_tot_compton = np.ma.concatenate(acInput_pair_flag_tot_compton)
	

	acInput_event_id_tot_compton = np.array(acInput_event_id_tot_compton)
	acInput_AC_panel_compton = np.array(acInput_AC_panel_compton)
	acInput_AC_subpanel_compton = np.array(acInput_AC_subpanel_compton)
	acInput_energy_dep_tot_compton = np.array(acInput_energy_dep_tot_compton)
	acInput_pair_flag_tot_compton = np.array(acInput_pair_flag_tot_compton)

	if (len(acInput_event_id_tot_ray) > 0):	
		acInput_event_id_tot_ray = np.ma.concatenate(acInput_event_id_tot_ray)
		acInput_AC_panel_ray = np.ma.concatenate(acInput_AC_panel_ray)
		acInput_AC_subpanel_ray = np.ma.concatenate(acInput_AC_subpanel_ray)
		acInput_energy_dep_tot_ray = np.ma.concatenate(acInput_energy_dep_tot_ray)
		acInput_pair_flag_tot_ray = np.ma.concatenate(acInput_pair_flag_tot_ray)	

	acInput_event_id_tot_ray = np.array(acInput_event_id_tot_ray)
	acInput_AC_panel_ray = np.array(acInput_AC_panel_ray)
	acInput_AC_subpanel_ray = np.array(acInput_AC_subpanel_ray)
	acInput_energy_dep_tot_ray = np.array(acInput_energy_dep_tot_ray)
	acInput_pair_flag_tot_ray = np.array(acInput_pair_flag_tot_ray)
	
	if (len(acInput_event_id_tot_pair) > 0):	
		acInput_event_id_tot_pair = np.ma.concatenate(acInput_event_id_tot_pair)
		acInput_AC_panel_pair = np.ma.concatenate(acInput_AC_panel_pair)
		acInput_AC_subpanel_pair = np.ma.concatenate(acInput_AC_subpanel_pair)
		acInput_energy_dep_tot_pair = np.ma.concatenate(acInput_energy_dep_tot_pair)
		acInput_pair_flag_tot_pair = np.ma.concatenate(acInput_pair_flag_tot_pair)
	

	acInput_event_id_tot_pair = np.array(acInput_event_id_tot_pair)
	acInput_AC_panel_pair = np.array(acInput_AC_panel_pair)
	acInput_AC_subpanel_pair = np.array(acInput_AC_subpanel_pair)
	acInput_energy_dep_tot_pair = np.array(acInput_energy_dep_tot_pair)
	acInput_pair_flag_tot_pair = np.array(acInput_pair_flag_tot_pair)
	
	
	S1_event_id_ac = np.ma.concatenate(S1_event_id_ac)
	S1_type_ac = np.ma.concatenate(S1_type_ac)
	S1_e_dep_ac = np.ma.concatenate(S1_e_dep_ac)
	S1_volume_id_ac = np.ma.concatenate(S1_volume_id_ac)
	S1_x_ac = np.ma.concatenate(S1_x_ac)
	S1_y_ac = np.ma.concatenate(S1_y_ac)
	S1_z_ac = np.ma.concatenate(S1_z_ac)
	
	S1_event_id_ac = np.array((S1_event_id_ac), dtype=np.int64)
	S1_type_ac = np.array((S1_type_ac), dtype=np.int64)
	S1_e_dep_ac = np.array((S1_e_dep_ac), dtype=np.float64)
	S1_volume_id_ac = np.array((S1_volume_id_ac), dtype=np.int64)
	S1_x_ac = np.array((S1_x_ac), dtype=np.float64)
	S1_y_ac = np.array((S1_y_ac), dtype=np.float64)
	S1_z_ac = np.array((S1_z_ac), dtype=np.float64)

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
		data.write('{:d}\t'.format(aa_strip_strip_id[r]))
		data.write('{:f}\t'.format(aa_strip_pos[r]))
		data.write('{:f}\t'.format(aa_strip_edep[r]))
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

	if os.path.exists(filepath+sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat'):
		os.remove(filepath+sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat')
		data = open(filepath+sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')
	else:
		data = open(filepath+sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')

	for r in range(len(aa_kalman_ray_event_id)):
		
		data.write('{:d}\t'.format(aa_kalman_ray_event_id[r]))
		data.write('{:d}\t'.format(aa_kalman_ray_theta_in[r]))
		data.write('{:d}\t'.format(aa_kalman_ray_phi_in[r]))
		data.write('{:s}\t'.format(aa_kalman_ray_ene_in[r]))
		data.write('{:d}\t'.format(aa_kalman_ray_plane_id[r]))
		data.write('{:f}\t'.format(aa_kalman_ray_zpos[r]))
		data.write('{:d}\t'.format(aa_kalman_ray_si_id[r]))
		data.write('{:f}\t'.format(aa_kalman_ray_pos[r]))
		data.write('{:f}\t'.format(aa_kalman_ray_edep[r]))
		data.write('{:d}\t'.format(aa_kalman_ray_strip_number[r]))
		data.write('{:d}\n'.format(aa_kalman_ray_pair[r]))

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
	col2 = fits.Column(name='VOL_ID', format='I', array=rawData_vol_id)
	col3 = fits.Column(name='MOTH_ID', format='J', array=rawData_moth_id)
	col4 = fits.Column(name='TRAY_ID', format='I', array=rawData_tray_id)
	col5 = fits.Column(name='PLANE_ID', format='I', array=rawData_plane_id)
	col6 = fits.Column(name='STRIP_ID_X', format='I', array=rawData_Strip_id_x)
	col7 = fits.Column(name='STRIP_ID_Y', format='I', array=rawData_Strip_id_y)
	col8 = fits.Column(name='E_DEP', format='F20.5', array=rawData_energy_dep)
	col9 = fits.Column(name='X_ENT', format='F20.5', array=rawData_ent_x)
	col10 = fits.Column(name='Y_ENT', format='F20.5', array=rawData_ent_y)
	col11 = fits.Column(name='Z_ENT', format='F20.5', array=rawData_ent_z)
	col12 = fits.Column(name='X_EXIT', format='F20.5', array=rawData_exit_x)
	col13 = fits.Column(name='Y_EXIT', format='F20.5', array=rawData_exit_y)
	col14 = fits.Column(name='Z_EXIT', format='F20.5', array=rawData_exit_z)
	col15 = fits.Column(name='PART_ID', format='I', array=rawData_part_id)
	col16 = fits.Column(name='TRK_ID', format='I', array=rawData_trk_id)
	col17 = fits.Column(name='CHILD_ID', format='I', array=rawData_child_id)
	col18 = fits.Column(name='PROC_ID', format='I', array=rawData_proc_id)
		
	cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18])
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
	col6 = fits.Column(name='TRK_FLAG', format='I', array=L0TRACKER_Glob_Si_id)
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
	col4 = fits.Column(name='TRK_FLAG', format='I', array=L05TRACKER_Glob_Si_id_cluster)
	col5 = fits.Column(name='POS', format='F20.5', array=L05TRACKER_Glob_pos_cluster)
	col6 = fits.Column(name='ZPOS', format='F20.5', array=L05TRACKER_Glob_zpos_cluster)
	col7 = fits.Column(name='E_DEP', format='F20.5', array=L05TRACKER_Glob_energy_cluster_dep)
	col8 = fits.Column(name='PAIR_FLAG', format='I', array=L05TRACKER_Glob_pair_flag_cluster)

		
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

	col1 = fits.Column(name='EVT_ID', format='IJ', array=rawData_event_id_cal)	
	col2 = fits.Column(name='E_DEP', format='1D', array=rawData_energy_dep_cal)
	col3 = fits.Column(name='X_ENT', format='1D', array=rawData_ent_x_cal)
	col4 = fits.Column(name='Y_ENT', format='1D', array=rawData_ent_y_cal)
	col5 = fits.Column(name='Z_ENT', format='1D', array=rawData_ent_z_cal)
	col6 = fits.Column(name='X_EXIT', format='1D', array=rawData_exit_x_cal)
	col7 = fits.Column(name='Y_EXIT', format='1D', array=rawData_exit_y_cal)
	col8 = fits.Column(name='Z_EXIT', format='1D', array=rawData_exit_z_cal)
	col9 = fits.Column(name='PART_ID', format='IJ', array=rawData_part_id_cal)
	col10 = fits.Column(name='TRK_ID', format='IJ', array=rawData_trk_id_cal)
	col11 = fits.Column(name='CHILD_ID', format='IJ', array=rawData_child_id_cal)
	col12 = fits.Column(name='PROC_ID', format='IJ', array=rawData_proc_id_cal)

		
	cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

	fits.setval(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy ='+ene_type, ext=1)

	fits.setval(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)

	fits.setval(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

	fits.setval(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Position unit = cm', ext=1)

	fits.setval(outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = keV', ext=1)



	col1 = fits.Column(name='EVT_ID', format='I', array=calInput_event_id_tot)	
	col2 = fits.Column(name='BAR_ID', format='I', array=calInput_bar_id_tot)
	col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=calInput_bar_ene_tot)
	col4 = fits.Column(name='PAIR_FLAG', format='I', array=calInput_pair_flag_tot)
		
	cols = fits.ColDefs([col1,col2,col3,col4])
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
	


	col1 = fits.Column(name='EVT_ID', format='I', array=calInput_event_id_tot_compton)	
	col2 = fits.Column(name='BAR_ID', format='I', array=calInput_bar_id_tot_compton)
	col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=calInput_bar_ene_tot_compton)
	col4 = fits.Column(name='PAIR_FLAG', format='I', array=calInput_pair_flag_tot_compton)
		
	cols = fits.ColDefs([col1,col2,col3,col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

	fits.setval(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

	fits.setval(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

	fits.setval(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

	fits.setval(outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)



	col1 = fits.Column(name='EVT_ID', format='I', array=calInput_event_id_tot_pair)	
	col2 = fits.Column(name='BAR_ID', format='I', array=calInput_bar_id_tot_pair)
	col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=calInput_bar_ene_tot_pair)
	col4 = fits.Column(name='PAIR_FLAG', format='I', array=calInput_pair_flag_tot_pair)
		
	cols = fits.ColDefs([col1,col2,col3,col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

	fits.setval(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

	fits.setval(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

	fits.setval(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

	fits.setval(outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)



	col1 = fits.Column(name='EVT_ID', format='I', array=calInput_event_id_tot_ray)	
	col2 = fits.Column(name='BAR_ID', format='I', array=calInput_bar_id_tot_ray)
	col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=calInput_bar_ene_tot_ray)
	col4 = fits.Column(name='PAIR_FLAG', format='I', array=calInput_pair_flag_tot_ray)
		
	cols = fits.ColDefs([col1,col2,col3,col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

	fits.setval(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

	fits.setval(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

	fits.setval(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

	fits.setval(outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)



	col1 = fits.Column(name='EVT_ID', format='I', array=calInputSum_event_id_tot)	
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

	col1 = fits.Column(name='EVT_ID', format='I', array=rawData_event_id_ac)	
	col2 = fits.Column(name='E_DEP', format='F20.5', array=rawData_energy_dep_ac)
	col3 = fits.Column(name='X_ENT', format='F20.5', array=rawData_ent_x_ac)
	col4 = fits.Column(name='Y_ENT', format='F20.5', array=rawData_ent_y_ac)
	col5 = fits.Column(name='Z_ENT', format='F20.5', array=rawData_ent_z_ac)
	col6 = fits.Column(name='X_EXIT', format='F20.5', array=rawData_exit_x_ac)
	col7 = fits.Column(name='Y_EXIT', format='F20.5', array=rawData_exit_y_ac)
	col8 = fits.Column(name='Z_EXIT', format='F20.5', array=rawData_exit_z_ac)
	col9 = fits.Column(name='PART_ID', format='I', array=rawData_part_id_ac)
	col10 = fits.Column(name='TRK_ID', format='I', array=rawData_trk_id_ac)
	col11 = fits.Column(name='CHILD_ID', format='I', array=rawData_child_id_ac)
	col12 = fits.Column(name='PROC_ID', format='I', array=rawData_proc_id_ac)
		
	cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')

	fits.setval(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

	fits.setval(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy ='+ene_type, ext=1)

	fits.setval(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)

	fits.setval(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

	fits.setval(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Position unit = cm', ext=1)

	fits.setval(outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = keV', ext=1)





	col1 = fits.Column(name='EVT_ID', format='I', array=acInput_event_id_tot)	
	col2 = fits.Column(name='AC_PANEL', format='A', array=acInput_AC_panel)
	col3 = fits.Column(name='AC_SUBPANEL', format='I', array=acInput_AC_subpanel)
	col4 = fits.Column(name='E_DEP', format='F20.15', array=acInput_energy_dep_tot)
	col4 = fits.Column(name='PAIR_FLAG', format='I', array=acInput_energy_dep_tot)

		
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
	


	col1 = fits.Column(name='EVT_ID', format='I', array=acInput_event_id_tot_compton)	
	col2 = fits.Column(name='AC_PANEL', format='A', array=acInput_AC_panel_compton)
	col3 = fits.Column(name='AC_SUBPANEL', format='I', array=acInput_AC_subpanel_compton)
	col4 = fits.Column(name='E_DEP', format='F20.15', array=acInput_energy_dep_tot_compton)
	col4 = fits.Column(name='PAIR_FLAG', format='I', array=acInput_energy_dep_tot_compton)

		
	cols = fits.ColDefs([col1,col2,col3,col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')


	fits.setval(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

	fits.setval(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

	fits.setval(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

	fits.setval(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

	fits.setval(filepath+'G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)



	col1 = fits.Column(name='EVT_ID', format='I', array=acInput_event_id_tot_pair)	
	col2 = fits.Column(name='AC_PANEL', format='A', array=acInput_AC_panel_pair)
	col3 = fits.Column(name='AC_SUBPANEL', format='I', array=acInput_AC_subpanel_pair)
	col4 = fits.Column(name='E_DEP', format='F20.15', array=acInput_energy_dep_tot_pair)
	col4 = fits.Column(name='PAIR_FLAG', format='I', array=acInput_energy_dep_tot_pair)

		
	cols = fits.ColDefs([col1,col2,col3,col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')


	fits.setval(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

	fits.setval(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

	fits.setval(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

	fits.setval(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

	fits.setval(filepath+'G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)


	col1 = fits.Column(name='EVT_ID', format='I', array=acInput_event_id_tot_ray)	
	col2 = fits.Column(name='AC_PANEL', format='A', array=acInput_AC_panel_ray)
	col3 = fits.Column(name='AC_SUBPANEL', format='I', array=acInput_AC_subpanel_ray)
	col4 = fits.Column(name='E_DEP', format='F20.15', array=acInput_energy_dep_tot_ray)
	col4 = fits.Column(name='PAIR_FLAG', format='I', array=acInput_energy_dep_tot_ray)

		
	cols = fits.ColDefs([col1,col2,col3,col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
	if os.path.exists(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits'):
		os.remove(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
		tbhdu.writeto(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')
	else:
		tbhdu.writeto(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits')


	fits.setval(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

	fits.setval(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

	fits.setval(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

	fits.setval(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

	fits.setval(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

	fits.setval(filepath+'G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+'all.fits', 'COMMENT', value='Energy unit = GeV', ext=1)



if isStrip == 1:

	if os.path.exists(filepath+sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat'):
		os.remove(filepath+sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat')

		data = open(filepath+sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')
	else:
		data = open(filepath+sim_tag+'_S1_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+'dat', 'w')


	if ((cal_flag == 0) & (ac_flag == 0)):

		for r in range(len(S1_event_id_trk)):
		
			data.write('{:10d}\t'.format(S1_event_id_trk[r]))
			data.write('{:10d}\t'.format(S1_type_trk[r]))
			data.write('{:10.5f}\t'.format(S1_e_dep_trk[r]))
			data.write('{:10d}\t'.format(S1_volume_id_trk[r]))
			data.write('{:10.5f}\t'.format(S1_x_trk[r]))
			data.write('{:10.5f}\t'.format(S1_y_trk[r]))
			data.write('{:10.5f}\n'.format(S1_z_trk[r]))
		
		data.close()

	
	if ((cal_flag == 1) & (ac_flag == 0)):
	
		S1_event_id = np.concatenate([S1_event_id_trk, S1_event_id_cal])
		S1_type = np.concatenate([S1_type_trk, S1_type_cal])
		S1_e_dep = np.concatenate([S1_e_dep_trk, S1_e_dep_cal])
		S1_volume_id = np.concatenate([S1_volume_id_trk, S1_volume_id_cal])
		S1_x = np.concatenate([S1_x_trk, S1_x_cal])
		S1_y = np.concatenate([S1_y_trk, S1_y_cal])
		S1_z = np.concatenate([S1_z_trk, S1_z_cal])
		
		id_sort = np.argsort(S1_event_id)
		S1_event_id_sort = S1_event_id[id_sort]
		S1_type_sort = S1_type[id_sort]
		S1_e_dep_sort = S1_e_dep[id_sort]
		S1_volume_id_sort = S1_volume_id[id_sort]
		S1_x_sort = S1_x[id_sort]
		S1_y_sort = S1_y[id_sort]
		S1_z_sort = S1_z[id_sort]
				
		for r in range(len(S1_event_id_sort)):
		
			data.write('{:10d}\t'.format(S1_event_id_sort[r]))
			data.write('{:10d}\t'.format(S1_type_sort[r]))
			data.write('{:10.5f}\t'.format(S1_e_dep_sort[r]))
			data.write('{:10d}\t'.format(S1_volume_id_sort[r]))
			data.write('{:10.5f}\t'.format(S1_x_sort[r]))
			data.write('{:10.5f}\t'.format(S1_y_sort[r]))
			data.write('{:10.5f}\n'.format(S1_z_sort[r]))
			
		data.close()

	if ((cal_flag == 0) & (ac_flag == 1)):
	
		S1_event_id = np.concatenate([S1_event_id_trk, S1_event_id_ac])
		S1_type = np.concatenate([S1_type_trk, S1_type_ac])
		S1_e_dep = np.concatenate([S1_e_dep_trk, S1_e_dep_ac])
		S1_volume_id = np.concatenate([S1_volume_id_trk, S1_volume_id_ac])
		S1_x = np.concatenate([S1_x_trk, S1_x_ac])
		S1_y = np.concatenate([S1_y_trk, S1_y_ac])
		S1_z = np.concatenate([S1_z_trk, S1_z_ac])
		
		id_sort = np.argsort(S1_event_id)
		S1_event_id_sort = S1_event_id[id_sort]
		S1_type_sort = S1_type[id_sort]
		S1_e_dep_sort = S1_e_dep[id_sort]
		S1_volume_id_sort = S1_volume_id[id_sort]
		S1_x_sort = S1_x[id_sort]
		S1_y_sort = S1_y[id_sort]
		S1_z_sort = S1_z[id_sort]
				
		for r in range(len(S1_event_id_sort)):
		
			data.write('{:10d}\t'.format(S1_event_id_sort[r]))
			data.write('{:10d}\t'.format(S1_type_sort[r]))
			data.write('{:10.5f}\t'.format(S1_e_dep_sort[r]))
			data.write('{:10d}\t'.format(S1_volume_id_sort[r]))
			data.write('{:10.5f}\t'.format(S1_x_sort[r]))
			data.write('{:10.5f}\t'.format(S1_y_sort[r]))
			data.write('{:10.5f}\n'.format(S1_z_sort[r]))
			
		data.close()

	
	if ((cal_flag == 1) & (ac_flag == 1)):


		S1_event_id = np.concatenate([S1_event_id_trk, S1_event_id_cal])
		S1_type = np.concatenate([S1_type_trk, S1_type_cal])
		S1_e_dep = np.concatenate([S1_e_dep_trk, S1_e_dep_cal])
		S1_volume_id = np.concatenate([S1_volume_id_trk, S1_volume_id_cal])
		S1_x = np.concatenate([S1_x_trk, S1_x_cal])
		S1_y = np.concatenate([S1_y_trk, S1_y_cal])
		S1_z = np.concatenate([S1_z_trk, S1_z_cal])
		
		S1_event_id = np.concatenate([S1_event_id, S1_event_id_ac])
		S1_type = np.concatenate([S1_type, S1_type_ac])
		S1_e_dep = np.concatenate([S1_e_dep, S1_e_dep_ac])
		S1_volume_id = np.concatenate([S1_volume_id, S1_volume_id_ac])
		S1_x = np.concatenate([S1_x, S1_x_ac])
		S1_y = np.concatenate([S1_y, S1_y_ac])
		S1_z = np.concatenate([S1_z, S1_z_ac])
		
		
		id_sort = np.argsort(S1_event_id)
		S1_event_id_sort = S1_event_id[id_sort]
		S1_type_sort = S1_type[id_sort]
		S1_e_dep_sort = S1_e_dep[id_sort]
		S1_volume_id_sort = S1_volume_id[id_sort]
		S1_x_sort = S1_x[id_sort]
		S1_y_sort = S1_y[id_sort]
		S1_z_sort = S1_z[id_sort]
				
		for r in range(len(S1_event_id_sort)):
		
			data.write('{:10d}\t'.format(S1_event_id_sort[r]))
			data.write('{:10d}\t'.format(S1_type_sort[r]))
			data.write('{:10.5f}\t'.format(S1_e_dep_sort[r]))
			data.write('{:10d}\t'.format(S1_volume_id_sort[r]))
			data.write('{:10.5f}\t'.format(S1_x_sort[r]))
			data.write('{:10.5f}\t'.format(S1_y_sort[r]))
			data.write('{:10.5f}\n'.format(S1_z_sort[r]))
			
		data.close()



