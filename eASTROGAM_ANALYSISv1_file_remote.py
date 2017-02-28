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
from function import *
from L05_cluster_baricenter import *
from summing_energy import *
from reading import *
from conversion import *
from writing import *
from Build_L0 import *
from L05_merging import *
from AC import *
from CAL import *
from flag_events import *
from Strip_analysis import *

sys.argv[0] = 'eASTROGAM_ANALYSISv1_file_remote.py'

##############################
#     parametri iniziali     #
##############################

astrogam_version = sys.argv[1]           # Enter eASTROGAM release (e.g. V1.0):
bogemms_tag = int(sys.argv[2])           # Enter BoGEMMS release (e.g. 211):
sim_type = int(sys.argv[3])              # Enter simulation type [0 = Mono, 1 = Range, 2 = Chen, 3: Vela, 4: Crab, 5: G400]:
py_list = int(sys.argv[4])               # Enter the Physics List [0 = QGSP_BERT_EMV, 100 = ARGO, 300 = FERMI, 400 = ASTROMEV]:
N_in = int(sys.argv[5])                  # Enter the number of emitted particles:
part_type = sys.argv[6]                  # Enter the particle type [ph = photons, mu = muons, g = geantino, p = proton, el = electron]:
ene_range = int(sys.argv[7])             # Enter energy distribution [0 = MONO, 1 = POW, 2 = EXP, 3 = LIN]:
ene_min = int(sys.argv[8])               # Enter miminum energy [MeV]:
ene_max = int(sys.argv[9])               # Enter maximum energy [MeV]:
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
ifile = int(sys.argv[22])		 # Enter the initial number of FITS files
n_fits = int(sys.argv[23])               # Enter the final number of FITS files :


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

#print(stripname,dir_cal)

# setting specific agile version variables

if astrogam_version=='V1.0':
	# --------> volume ID
	tracker_top_vol_start = 1090000
	tracker_bottom_vol_start = 1000000
	tracker_top_bot_diff = 90000

	cal_vol_start = 50000
	cal_vol_end = 58463

	ac_vol_start = 301
	ac_vol_end = 350

	panel_S = [301, 302, 303]
	panel_D = [311, 312, 313]
	panel_F = [321, 322, 323]
	panel_B = [331, 332, 333]
	panel_top = 340

	# --------> design
	N_tray = 56
	N_plane = N_tray*1
	N_strip = 3840
	tray_side = 92.16 #cm
	strip_side = tray_side/N_strip

	# --------> processing
	# accoppiamento capacitivo
	# acap = [0.035, 0.045, 0.095, 0.115, 0.38, 1., 0.38, 0.115, 0.095, 0.045, 0.035]
	# tracker energy threshold (0.25 MIP)
	E_th = float(energy_thresh)  # keV

	E_th_cal = 30. # keV


#parte lettura file fits


filepath = './input_eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+dir_cal+dir_passive+'/'+ene_type+'MeV.'+sim_name+'.'+str(theta_type)+'theta.'+pol_string+str(N_in)+part_type

print('eASTROGAM simulation path:' + filepath)

outdir = ('./output_eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+str(sim_name)+'/'+str(ene_type)+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV')

if not os.path.exists(outdir):
	out_dir = os.makedirs(outdir,0777)

#ifile = 0
while ifile <= n_fits:
	start = time.time()
	
	print('Reading the THELSIM file.....'+ str(ifile))

	vol_id_tr, moth_id_tr, event_id_tr, en_dep_tr, x_en_tr, y_en_tr, z_en_tr, x_ex_tr, y_ex_tr, z_ex_tr, child_id_tr, proc_id_tr, theta_ent_tr, phi_ent_tr, theta_exit_tr, phi_exit_tr, x_pos, y_pos, z_pos, vol_id_cal, moth_id_cal, event_id_cal, energy_dep_cal, x_en_cal, y_en_cal, z_en_cal, x_ex_cal, y_ex_cal, z_ex_cal, child_id_cal, proc_id_cal, theta_ent_cal, phi_ent_cal, theta_exit_cal, phi_exit_cal, vol_id_ac, event_id_ac, en_dep_ac, x_en_ac, y_en_ac, z_en_ac, x_ex_ac, y_ex_ac, z_ex_ac, child_id_ac, proc_id_ac, theta_ent_ac, phi_ent_ac, theta_exit_ac, phi_exit_ac = reading_fits(filepath, ifile, cal_flag, ac_flag, tracker_top_vol_start, tracker_bottom_vol_start, tracker_top_bot_diff, cal_vol_start, cal_vol_end, ac_vol_start, ac_vol_end, part_type)
	
	
	
			#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			#%                             Processing the tracker                          %
			#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if astrogam_version == 'V1.0':
		# From Tracker volume ID to strip and tray ID and conversion from tray ID (starting from bottom) to plane ID (starting from the top)

		Strip_id_x, Strip_id_y, tray_id, plane_id = conversion(vol_id_tr, isStrip, repli, moth_id_tr, tracker_bottom_vol_start, N_tray, tracker_top_bot_diff)	

		tray_id = np.array(tray_id)
		plane_id = np.array(plane_id)
		Strip_id_x = np.array(Strip_id_x)
		Strip_id_y = np.array(Strip_id_y)


	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('                             Tracker   '                     )
	print('           Saving the Tracker raw hits (fits and .dat)      ')
	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	if astrogam_version == 'V1.0':
		
		writing_G4raw(event_id_tr, vol_id_tr, moth_id_tr, tray_id, plane_id, Strip_id_x, Strip_id_y, en_dep_tr, x_en_tr, y_en_tr, z_en_tr, x_ex_tr, y_ex_tr, z_ex_tr, child_id_tr, proc_id_tr, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile)		


		if isStrip == 0:
	
			# ASCII Columns:
			# - c1 = event ID
			# - c2 = theta input
			# - c3 = phi input
			# - c4 = energy input
			# - c5 = plane ID
			# - c6 = Pos Z
			# - c7 = X/Y flag (X = 0, Y = 1)
			# - c8 = Cluster position (reference system center at the Silicon layer center)
			# - c9 = energy deposition (keV)
			# - c10 = number of strips composing the cluster
			# - c11 = child id
			# - c12 = proc id 
			
			writing_AA_fake(event_id_tr, plane_id, x_pos, y_pos, z_pos, en_dep_tr, child_id_tr, proc_id_tr, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile)		
				
		
		# Loading the LUT

		if isStrip == 1:
			
			filename_x_top = './conf/ARCH.XSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS'
			filename_y_top = './conf/ARCH.YSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS'

			struct_x_top = fits.open(filename_x_top)
			struct_y_top = fits.open(filename_y_top)

			g_x_data = struct_x_top[1].data
			
			Arch_vol_id_x_top = g_x_data.field('VOLUME_ID')
			Arch_moth_id_x_top = g_x_data.field('MOTHER_ID')
			Arch_Strip_id_x_top = g_x_data.field('STRIP_ID')
			Arch_Si_id_x_top = g_x_data.field('TRK_FLAG')
			Arch_tray_id_x_top = g_x_data.field('TRAY_ID')
			Arch_plane_id_x_top = g_x_data.field('PLANE_ID')
			Arch_xpos_x_top = g_x_data.field('XPOS')
			Arch_zpos_x_top = g_x_data.field('ZPOS')
			Arch_energy_dep_x_top = g_x_data.field('E_DEP')
			Arch_pair_flag_x_top = g_x_data.field('PAIR_FLAG')

			g_y_data = struct_y_top[1].data

			Arch_vol_id_y_top = g_y_data.field('VOLUME_ID')
			Arch_moth_id_y_top = g_y_data.field('MOTHER_ID')
			Arch_Strip_id_y_top = g_y_data.field('STRIP_ID')
			Arch_Si_id_y_top = g_y_data.field('TRK_FLAG')
			Arch_tray_id_y_top = g_y_data.field('TRAY_ID')
			Arch_plane_id_y_top = g_y_data.field('PLANE_ID')
			Arch_ypos_y_top = g_y_data.field('YPOS')
			Arch_zpos_y_top = g_y_data.field('ZPOS')
			Arch_energy_dep_y_top = g_y_data.field('E_DEP')
			Arch_pair_flag_y_top = g_y_data.field('PAIR_FLAG')
		

			e_dep_temp, event_id_tot, vol_id_tot, moth_id_tot, Strip_id_x_tot, Strip_id_y_tot, tray_id_tot, plane_id_tot, energy_dep_tot, pair_flag_tot = flag_events(event_id_tr, vol_id_tr, moth_id_tr, Strip_id_x, Strip_id_y, tray_id, plane_id, en_dep_tr, child_id_tr, proc_id_tr)


			event_id_tot_temp = np.zeros(2*len(event_id_tot), dtype = np.int64)
			vol_id_tot_temp = np.zeros(2*len(event_id_tot), dtype = np.int64)
			moth_id_tot_temp = np.zeros(2*len(event_id_tot), dtype = np.int64)
			Strip_id_tot_temp = np.zeros(2*len(event_id_tot), dtype = np.int64)
			Si_id_tot_temp = np.zeros(2*len(event_id_tot), dtype = np.int64)
			tray_id_tot_temp = np.zeros(2*len(event_id_tot), dtype = np.int64)
			plane_id_tot_temp = np.zeros(2*len(event_id_tot), dtype = np.int64)
			energy_dep_tot_temp = np.zeros(2*len(event_id_tot))
			pair_flag_tot_temp = np.zeros(2*len(event_id_tot))

			jev = 0
			while jev < len(event_id_tot):
				ev_index = jev * 2

				
				event_id_tot_temp[ev_index] = event_id_tot[jev]
				event_id_tot_temp[ev_index+1] = event_id_tot[jev]
				vol_id_tot_temp[ev_index] = Strip_id_x_tot[jev]
				vol_id_tot_temp[ev_index+1] = Strip_id_y_tot[jev]
				moth_id_tot_temp[ev_index] = moth_id_tot[jev] - Strip_id_x_tot[jev]
				moth_id_tot_temp[ev_index+1] = moth_id_tot[jev] - Strip_id_x_tot[jev]
				Strip_id_tot_temp[ev_index] = Strip_id_x_tot[jev]
				Strip_id_tot_temp[ev_index+1] = Strip_id_y_tot[jev]
				Si_id_tot_temp[ev_index] = 0
				Si_id_tot_temp[ev_index+1] = 1
				tray_id_tot_temp[ev_index] = tray_id_tot[jev]
				tray_id_tot_temp[ev_index+1] = tray_id_tot[jev]
				plane_id_tot_temp[ev_index] = plane_id_tot[jev]
				plane_id_tot_temp[ev_index+1] = plane_id_tot[jev]
				energy_dep_tot_temp[ev_index] = energy_dep_tot[jev]/2.
				energy_dep_tot_temp[ev_index+1] = energy_dep_tot[jev]/2.
				pair_flag_tot_temp[ev_index] = pair_flag_tot[jev]
				pair_flag_tot_temp[ev_index+1] = pair_flag_tot[jev]

				
				jev = jev + 1
			
			
			

			#
			# Summing the energy along the strip and applying the energy thresold
			#
				
			event_id_tot, vol_id_tot, moth_id_tot, Strip_id_tot, Si_id_tot, tray_id_tot, plane_id_tot, energy_dep_tot, pair_flag_tot = summing_energy(E_th, event_id_tot_temp, vol_id_tot_temp, moth_id_tot_temp, Strip_id_tot_temp, Si_id_tot_temp, tray_id_tot_temp, plane_id_tot_temp, energy_dep_tot_temp, pair_flag_tot_temp)

			####Funzione indici unici event_id_tot

			event_uniq_index = index_uniq(event_id_tot)

			#### Fine Funzione indici unici event_id_tot
		
			N_trig = len(event_uniq_index)

			event_array = event_id_tot[event_uniq_index]


	if astrogam_version == 'V1.0':

		if isStrip == 1:

			#Total number of strips
			Total_vol_x_top = (N_tray)*N_strip
			Total_vol_y_top = (N_tray)*N_strip


			print('Number of tracker triggered events: '+ str(N_trig))
			

			Glob_vol_id_x_top, Glob_moth_id_x_top, Glob_Strip_id_x_top, Glob_Si_id_x_top, Glob_tray_id_x_top, Glob_plane_id_x_top, Glob_xpos_x_top, Glob_zpos_x_top, Glob_energy_dep_x_top, Glob_pair_flag_x_top, Glob_vol_id_y_top, Glob_moth_id_y_top, Glob_Strip_id_y_top, Glob_Si_id_y_top, Glob_tray_id_y_top, Glob_plane_id_y_top, Glob_ypos_y_top, Glob_zpos_y_top, Glob_energy_dep_y_top, Glob_pair_flag_y_top, N_ev = strip_analysis(Total_vol_x_top, Total_vol_y_top, N_trig, Arch_vol_id_x_top, Arch_moth_id_x_top, Arch_Strip_id_x_top, Arch_Si_id_x_top, Arch_tray_id_x_top, Arch_plane_id_x_top, Arch_xpos_x_top, Arch_zpos_x_top, Arch_energy_dep_x_top, Arch_pair_flag_x_top, Arch_vol_id_y_top, Arch_moth_id_y_top, Arch_Strip_id_y_top, Arch_Si_id_y_top, Arch_tray_id_y_top, Arch_plane_id_y_top, Arch_ypos_y_top, Arch_zpos_y_top, Arch_energy_dep_y_top, Arch_pair_flag_y_top, event_id_tot, vol_id_tot, moth_id_tot, Strip_id_tot, Si_id_tot, tray_id_tot, plane_id_tot, energy_dep_tot, pair_flag_tot)
		
			print('N_ev: '+ str(N_ev))				


			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('              Build the LEVEL 0 output            ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')		

			Glob_event_id_test, Glob_vol_id_test, Glob_moth_id_test, Glob_Strip_id_test, Glob_Si_id_test, Glob_tray_id_test, Glob_plane_id_test, Glob_pos_test, Glob_zpos_test, Glob_energy_dep_test, Glob_pair_flag_test = build_L0(Glob_energy_dep_x_top, Glob_energy_dep_y_top, Glob_vol_id_x_top, Glob_vol_id_y_top, Glob_moth_id_x_top, Glob_moth_id_y_top, Glob_Strip_id_x_top, Glob_Strip_id_y_top, Glob_Si_id_x_top, Glob_Si_id_y_top, Glob_tray_id_x_top, Glob_tray_id_y_top, Glob_plane_id_x_top, Glob_plane_id_y_top, Glob_xpos_x_top, Glob_ypos_y_top, Glob_zpos_x_top, Glob_zpos_y_top, Glob_pair_flag_x_top, Glob_pair_flag_y_top, N_trig, event_array)

			# Level 0 = energy summed
			# Level 0 = the events are sorted in tray, and Y before X within the same tray
			# energy threshold applied		
			
			writing_L0(Glob_event_id_test, Glob_vol_id_test, Glob_moth_id_test, Glob_tray_id_test, Glob_plane_id_test, Glob_Si_id_test, Glob_Strip_id_test, Glob_pos_test, Glob_zpos_test, Glob_energy_dep_test, Glob_pair_flag_test, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile, N_trig) 



			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - strip')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

			# ASCII Columns:
			# - c1 = event ID
			# - c2 = theta input
			# - c3 = phi input
			# - c4 = energy input
			# - c5 = plane ID
			# - c6 = Pos Z
			# - c7 = X/Y flag (X = 0, Y = 1)
			# - c8 = Strip ID
			# - c9 = strip position (reference system center at the Silicon layer center)
			# - c10 = energy deposition (keV)
			# - c11 = pair flag (1 = pair, 0 = not pair)

			writing_AA_strip(Glob_event_id_test, Glob_Si_id_test, Glob_tray_id_test, Glob_plane_id_test, Glob_Strip_id_test, Glob_pos_test, Glob_zpos_test, Glob_energy_dep_test, Glob_pair_flag_test, outdir, sim_tag, N_in, part_type, sname, ene_dis, ang_type, ene_type, theta_type, phi_type, pol_string, ifile)


			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('       L0.5 - cluster baricenter ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


			print('N_trig: '+ str(N_trig))


#			Glob_event_id_x_top_cluster = []
#			Glob_Si_id_x_top_cluster = []
#			Glob_tray_id_x_top_cluster = []
#			Glob_plane_id_x_top_cluster = []
#			Glob_zpos_x_top_cluster = []
#			Glob_energy_dep_x_top_cluster = []
#			Glob_xpos_x_top_cluster = []
#			Glob_Strip_number_x_top_cluster = []
#			Glob_pair_flag_x_top_cluster = []
			
			Glob_event_id_x_top_cluster, Glob_Si_id_x_top_cluster, Glob_tray_id_x_top_cluster, Glob_plane_id_x_top_cluster, Glob_zpos_x_top_cluster, Glob_energy_dep_x_top_cluster, Glob_xpos_x_top_cluster, Glob_Strip_number_x_top_cluster, Glob_pair_flag_x_top_cluster = L05_cluster_x(N_trig, Glob_plane_id_x_top, Glob_vol_id_x_top, Glob_moth_id_x_top, Glob_Strip_id_x_top, Glob_Si_id_x_top, Glob_tray_id_x_top, Glob_xpos_x_top, Glob_zpos_x_top, Glob_energy_dep_x_top, Glob_pair_flag_x_top)


			#Glob_event_id_y_top_cluster = []
			#Glob_Si_id_y_top_cluster = []
			#Glob_tray_id_y_top_cluster = []
			#Glob_plane_id_y_top_cluster = []
			#Glob_zpos_y_top_cluster = []
			#Glob_energy_dep_y_top_cluster = []
			#Glob_ypos_y_top_cluster = []
			#Glob_Strip_number_y_top_cluster = []
			#Glob_pair_flag_y_top_cluster = []

			Glob_event_id_y_top_cluster, Glob_Si_id_y_top_cluster, Glob_tray_id_y_top_cluster, Glob_plane_id_y_top_cluster, Glob_zpos_y_top_cluster, Glob_energy_dep_y_top_cluster, Glob_ypos_y_top_cluster, Glob_Strip_number_y_top_cluster, Glob_pair_flag_y_top_cluster = L05_cluster_y(N_trig, Glob_plane_id_y_top, Glob_vol_id_y_top, Glob_moth_id_y_top, Glob_Strip_id_y_top, Glob_Si_id_y_top, Glob_tray_id_y_top, Glob_ypos_y_top, Glob_zpos_y_top, Glob_energy_dep_y_top, Glob_pair_flag_y_top)


			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('             L0.5 - X-Y layers merging ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

			Glob_event_id_cluster, Glob_Si_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_Strip_number_cluster, Glob_pair_flag_cluster = merging(event_array, N_trig, Glob_event_id_x_top_cluster, Glob_event_id_y_top_cluster, Glob_Strip_number_x_top_cluster, Glob_Strip_number_y_top_cluster, Glob_Si_id_x_top_cluster, Glob_Si_id_y_top_cluster, Glob_tray_id_x_top_cluster, Glob_tray_id_y_top_cluster, Glob_plane_id_x_top_cluster, Glob_plane_id_y_top_cluster, Glob_xpos_x_top_cluster, Glob_ypos_y_top_cluster, Glob_zpos_x_top_cluster, Glob_zpos_y_top_cluster, Glob_energy_dep_x_top_cluster, Glob_energy_dep_y_top_cluster, Glob_pair_flag_x_top_cluster, Glob_pair_flag_y_top_cluster)			

			# Level 0.5 = energy summed, MIP threshold applied, strip position used

			writing_L05(Glob_event_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_Si_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_pair_flag_cluster, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile, N_trig)


			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - cluster')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			
			# ASCII Columns:
			# - c1 = event ID
			# - c2 = theta input
			# - c3 = phi input
			# - c4 = energy input
			# - c5 = plane ID
			# - c6 = Pos Z
			# - c7 = X/Y flag (X = 0, Y = 1)
			# - c8 = Cluster position (reference system center at the Silicon layer center)
			# - c9 = energy deposition (keV)
			# - c10 = number of strips composing the cluster
			# - c11 = pair flag (1 = pair, 0 = not pair, 2 = compton)

			writing_AA_cluster(Glob_event_id_cluster, Glob_Si_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_pair_flag_cluster, Glob_Strip_number_cluster, outdir, sim_tag, N_in, part_type, sname, ene_dis, ang_type, ene_type, theta_type, phi_type, pol_string, ifile)
					
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - pairs')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			
			# ASCII Columns:
			# - c1 = event ID
			# - c2 = theta input
			# - c3 = phi input
			# - c4 = energy input
			# - c5 = plane ID
			# - c6 = Pos Z
			# - c7 = X/Y flag (X = 0, Y = 1)
			# - c8 = Cluster position (reference system center at the Silicon layer center)
			# - c9 = energy deposition (keV)
			# - c10 = number of strips composing the cluster
			# - c11 = pair flag (1 = pair)

			writing_AA_cluster_pairs(Glob_event_id_cluster, Glob_Si_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_pair_flag_cluster, Glob_Strip_number_cluster, outdir, sim_tag, N_in, part_type, sname, ene_dis, ang_type, ene_type, theta_type, phi_type, pol_string, ifile)

			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - compton')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
											   
			# ASCII Columns:
			# - c1 = event ID
			# - c2 = theta input
			# - c3 = phi input
			# - c4 = energy input
			# - c5 = plane ID
			# - c6 = Pos Z
			# - c7 = X/Y flag (X = 0, Y = 1)
			# - c8 = Cluster position (reference system center at the Silicon layer center)
			# - c9 = energy deposition (keV)
			# - c10 = number of strips composing the cluster
			# - c11 = pair flag (2 = compton)
									   
			writing_AA_cluster_compton(Glob_event_id_cluster, Glob_Si_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_pair_flag_cluster, Glob_Strip_number_cluster, outdir, sim_tag, N_in, part_type, sname, ene_dis, ang_type, ene_type, theta_type, phi_type, pol_string, ifile)



	if cal_flag == 1:
		

		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		#                             Processing the calorimeter
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('          Calorimeter Bar Energy attenuation'        )                
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			
		bar_ene = energy_dep_cal

		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('                   Calorimeter                ')
		print('              Applying the minimum cut                ')
		print('                Summing the energy                ')
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

		event_id_tot_cal, bar_id_tot, bar_ene_tot = G4_cal(event_id_cal, vol_id_cal, moth_id_cal, bar_ene, E_th_cal, cal_vol_start)
		
		writing_G4cal(event_id_tot_cal, bar_id_tot, bar_ene_tot, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile)

		event_id_tot_cal_sum, bar_ene_tot_sum = cal_sum(event_id_tot_cal, bar_ene_tot)
		
		writing_cal_sum(event_id_tot_cal_sum, bar_ene_tot_sum, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile)
		



	if ac_flag == 1:

		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('                          AC')
		print('                  Summing the energy                ')
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

		AC_panel, AC_subpanel, event_id_tot_ac, energy_dep_tot_ac = AC_analysis(event_id_ac, vol_id_ac, energy_dep_ac, panel_S, panel_D, panel_F, panel_B, panel_top)

		writing_G4ac(event_id_tot_ac, AC_panel, AC_subpanel, energy_dep_tot_ac, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile)

############    passo a file fits successivo   ####################

	end = time.time() - start

	print('Time: ' + str(end))
				
	ifile = ifile + 1



	













