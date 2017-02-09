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

sys.argv[0] = 'eASTROGAM_ANALYSISv1_file_remote.py'
##############################
#         Funzioni           #
##############################

########## index_uniq da utilizzare solo con array gia' ordinati ########
def index_uniq(a):							#
									#
	a_index = []							#
									#
	for b in range(len(a)):						#
		a_index_old = b						#
		a_index.append(a_index_old)				#
									#
	a_uniq_index = []						#
									#
	b = 1								#
	while b < len(a):						#
									#
		if a[b] != a[b-1]:					#
			a_uniq_index.append(a_index[b-1])		#
									#
		if b == len(a)-1:					#
			a_uniq_index.append(a_index[b])			#
									#
		b = b + 1						#
									#
	return(a_uniq_index)						#
									#
#########################################################################



##############################
#     parametri iniziali     #
##############################

astrogam_version = sys.argv[1]           # Enter eASTROGAM release (e.g. V1.0):
bogemms_tag = int(sys.argv[2])           # Enter BoGEMMS release (e.g. 211):
sim_type = int(sys.argv[3])              # Enter simulation type [0 = Mono, 1 = Range, 2 = Chen, 3: Vela, 4: Crab, 5: G400]:
py_list = int(sys.argv[4])               # Enter the Physics List [0 = QGSP_BERT_EMV, 100 = ARGO, 300 = FERMI, 400 = ASTROMEV]:
N_in = int(sys.argv[5])                  # Enter the number of emitted particles:
part_type = sys.argv[6]                  # Enter the particle type [ph = photons, mu = muons, g = geantino, p = proton, el = electron]:
n_fits = int(sys.argv[7])                # Enter number of FITS files:
ene_range = int(sys.argv[8])             # Enter energy distribution [0 = MONO, 1 = POW, 2 = EXP, 3 = LIN]:
ene_min = int(sys.argv[9])               # Enter miminum energy [MeV]:
ene_max = int(sys.argv[10])              # Enter maximum energy [MeV]:
ang_type = sys.argv[11]                  # Enter the angular distribution [e.g. UNI, ISO]:
theta_type = int(sys.argv[12])           # Enter theta:
phi_type = int(sys.argv[13])             # Enter phi:
pol_type = int(sys.argv[14])             # Is the source polarized? [0 = false, 1 = true]:
pol_angle = int(sys.argv[15])            # Enter Polarization angle:
source_g = int(sys.argv[16])             # Enter source geometry [0 = Point, 1 = Plane]:
isStrip = int(sys.argv[17])              # Strip/Pixels activated? [0 = false, 1 = true]
repli = int(sys.argv[18])                # Strips/Pixels replicated? [0 = false, 1 = true]
cal_flag = int(sys.argv[19])             # Is Cal present? [0 = false, 1 = true]:
ac_flag = int(sys.argv[20])              # Is AC present? [0 = false, 1 = true]:
passive_flag = int(sys.argv[21])         # Is Passive present? [0 = false, 1 = true]:
energy_thresh = int(sys.argv[22])        # Enter energy threshold [keV]:
ifile = int(sys.argv[23])		 # Enter the initial file number


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
	N_tray = 56l
	N_plane = N_tray*1
	N_strip = 3840l
	tray_side = 92.16 #cm
	strip_side = tray_side/N_strip

	# --------> processing
	# accoppiamento capacitivo
	# acap = [0.035, 0.045, 0.095, 0.115, 0.38, 1., 0.38, 0.115, 0.095, 0.045, 0.035]
	# tracker energy threshold (0.25 MIP)
	E_th = float(energy_thresh)  # keV

	E_th_cal = 30. # keV


#parte lettura file fits


filepath = '/eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+dir_cal+dir_passive+'/'+ene_type+'MeV.'+sim_name+'.'+str(theta_type)+'theta.'+pol_string+str(N_in)+part_type

print('eASTROGAM simulation path:' + filepath)

outdir = ('./eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+str(sim_name)+'/'+str(ene_type)+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV')

if not os.path.exists('./eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+str(sim_name)+'/'+str(ene_type)+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV'):
	out_dir = os.makedirs(outdir,0777)

#ifile = 0
while ifile < n_fits:
	start = time.time()
	
	print('Reading the THELSIM file.....'+ str(ifile))

	t = fits.open('./eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+dir_cal+dir_passive+'/'+ene_type+'MeV.'+sim_name+'.'+str(theta_type)+'theta.'+pol_string+str(N_in)+part_type+'/xyz.'+str(ifile)+'.fits.gz')
   
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
				

		# Reading the Calorimeter (controllare separazione geantino e altro)
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

	
	
	
			#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			#%                             Processing the tracker                          %
			#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if astrogam_version == 'V1.0':
		# From Tracker volume ID to strip and tray ID	
		Strip_id_x = []
		Strip_id_y = []
		tray_id = []
		
		# Conversion from tray ID (starting from bottom) to plane ID (starting from the top)
		plane_id = []

		j=0
		
		while j < len(vol_id_tr):
		
			if isStrip == 1:
				if repli == 1:
					Strip_y = vol_id_tr[j]
					tray = moth_id_tr[j]/tracker_bottom_vol_start
					invert_tray_id = (N_tray - tray)+1
					vol_id_temp = moth_id_tr[j] - (tracker_bottom_vol_start*tray + tracker_top_bot_diff) # removing 1000000xn_tray + 90000					#
					Strip_x = vol_id_temp
					plane = invert_tray_id	

					plane_id.append(plane)
					Strip_id_y.append(Strip_y)
					Strip_id_x.append(Strip_x)
					tray_id.append(tray)

			else:	
					

				Strip_y = 0
				tray = vol_id_tr[j]/tracker_bottom_vol_start
				invert_tray_id = (N_tray - tray)+1
				Strip_x= 0
				plane = invert_tray_id				
				
				plane_id.append(plane)
				Strip_id_y.append(Strip_y)
				Strip_id_x.append(Strip_x)
				tray_id.append(tray)

			j = j+1


	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('                             Tracker   '                     )
	print('           Saving the Tracker raw hits (fits and .dat)      ')
	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	if astrogam_version == 'V1.0':
 	
		event_id_tr = np.array(event_id_tr)
		vol_id_tr = np.array(vol_id_tr)
		moth_id_tr = np.array(moth_id_tr)
		tr_id = np.array(tray_id)
		pl_id = np.array(plane_id)
		Strip_id_x = np.array(Strip_id_x)
		Strip_id_y = np.array(Strip_id_y)
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

	
		col1 = fits.Column(name='EVT_ID', format='I', array=event_id_tr)	
		col2 = fits.Column(name='VOL_ID', format='I', array=vol_id_tr)
		col3 = fits.Column(name='MOTH_ID', format='J', array=moth_id_tr)
		col4 = fits.Column(name='TRAY_ID', format='I', array=tr_id)
		col5 = fits.Column(name='PLANE_ID', format='I', array=pl_id)
		col6 = fits.Column(name='STRIP_ID_X', format='I', array=Strip_id_x)
		col7 = fits.Column(name='STRIP_ID_Y', format='I', array=Strip_id_y)
		col8 = fits.Column(name='E_DEP', format='F20.5', array=en_dep_tr)
		col9 = fits.Column(name='X_ENT', format='F20.5', array=x_en_tr)
		col10 = fits.Column(name='Y_ENT', format='F20.5', array=y_en_tr)
		col11 = fits.Column(name='Z_ENT', format='F20.5', array=z_en_tr)
		col12 = fits.Column(name='X_EXIT', format='F20.5', array=x_ex_tr)
		col13 = fits.Column(name='Y_EXIT', format='F20.5', array=y_en_tr)
		col14 = fits.Column(name='Z_EXIT', format='F20.5', array=z_en_tr)
		col15 = fits.Column(name='CHILD_ID', format='I', array=child_id_tr)
		col16 = fits.Column(name='PROC_ID', format='I', array=proc_id_tr)
		
		cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+str(ifile)+'.fits'):
			os.remove(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+str(ifile)+'.fits')
			tbhdu.writeto(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')

		fits.setval(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

		fits.setval(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy ='+ene_type, ext=1)

		fits.setval(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)

		fits.setval(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

		fits.setval(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Position unit = cm', ext=1)

		fits.setval(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = keV', ext=1)


		if isStrip == 0:
			if os.path.exists(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
				os.remove(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
				data = open(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
			else:
				data = open(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
	
			

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


			j = 0
			while 1:
				where_event_eq = np.where(event_id_tr == event_id_tr[j])
				where_event_eq = where_event_eq[0]

				plane_id_temp = pl_id[where_event_eq]
				Cluster_x_temp  = x_pos[where_event_eq]
				Cluster_y_temp  = y_pos[where_event_eq]
				Cluster_z_temp  = z_pos[where_event_eq]
				e_dep_x_temp  = (en_dep_tr[where_event_eq])/2.
				e_dep_y_temp  = (en_dep_tr[where_event_eq])/2.
				child_temp = child_id_tr[where_event_eq]
				proc_temp = proc_id_tr[where_event_eq]

				# ------------------------------------
				
				# X VIEW

				r = 0
				while r < len(Cluster_x_temp):

					if e_dep_x_temp[r] > E_th:
						data.write('{:d}\t'.format(event_id_tr[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(ene_type))
						data.write('{:d}\t'.format(plane_id_temp[r]))
						data.write('{:f}\t'.format(Cluster_z_temp[r]))
						data.write('{:d}\t'.format(0))
						data.write('{:f}\t'.format(Cluster_x_temp[r]))
						data.write('{:f}\t'.format(e_dep_x_temp[r]))
						data.write('{:d}\t'.format(1))
						data.write('{:d}\t'.format(child_temp[r]))
						data.write('{:d}\n'.format(proc_temp[r]))

					r = r + 1
				# ------------------------------------

				# Y VIEW

				r = 0
				while r < len(Cluster_y_temp):
					if e_dep_y_temp[r] > E_th:
						data.write('{:d}\t'.format(event_id_tr[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(ene_type))
						data.write('{:d}\t'.format(plane_id_temp[r]))
						data.write('{:f}\t'.format(Cluster_z_temp[r]))
						data.write('{:d}\t'.format(1))
						data.write('{:f}\t'.format(Cluster_y_temp[r]))
						data.write('{:f}\t'.format(e_dep_y_temp[r]))
						data.write('{:d}\t'.format(1))
						data.write('{:d}\t'.format(child_temp[r]))
						data.write('{:d}\n'.format(proc_temp[r]))

					r = r + 1


				N_event_eq = len(where_event_eq)                                          
				if where_event_eq[N_event_eq-1] < len(event_id_tr)-1:
					j = where_event_eq[N_event_eq-1]+1					
				else:
					break

			data.close()			
				
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
				tray_id_temp  = tr_id[where_event_eq]
				plane_id_temp  = pl_id[where_event_eq]
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
			# Summing the energy along the strip
			#

			e_dep_temp = [] 		
			event_id_tot = []
			vol_id_tot = []
			moth_id_tot = []
			Strip_id_tot = []
			Si_id_tot = []
			tray_id_tot = []
			plane_id_tot = []
			energy_dep_tot = []
			pair_flag_tot = []


			j = 0
			while 1:           
				where_event_eq = np.where(event_id_tot_temp == event_id_tot_temp[j])
				where_event_eq = where_event_eq[0]


				vol_id_temp = vol_id_tot_temp[where_event_eq]
				moth_id_temp = moth_id_tot_temp[where_event_eq]
				Strip_id_temp = Strip_id_tot_temp[where_event_eq]
				Si_id_temp = Si_id_tot_temp[where_event_eq]
				tray_id_temp = tray_id_tot_temp[where_event_eq]
				plane_id_temp = plane_id_tot_temp[where_event_eq]
				energy_dep_temp = energy_dep_tot_temp[where_event_eq]
				pair_flag_temp = pair_flag_tot_temp[where_event_eq]
				

				
				r = 0												
				while 1:

					where_vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]) & (Si_id_temp == 0))	
					where_vol_eq = where_vol_eq[0]
					
					where_other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]) | (Si_id_temp != 0))
					where_other_vol = where_other_vol[0]

					if len(where_vol_eq) != 0:
						e_dep_temp_old = np.sum(energy_dep_temp[where_vol_eq])
						event_id_tot_old = event_id_tot_temp[j]
						vol_id_tot_old = vol_id_temp[r]
						moth_id_tot_old = moth_id_temp[r]
						Si_id_tot_old = 0
						Strip_id_tot_old = Strip_id_temp[r]
						tray_id_tot_old = tray_id_temp[r]
						plane_id_tot_old = plane_id_temp[r]
						energy_dep_tot_old = e_dep_temp_old
						pair_flag_tot_old = pair_flag_temp[r]

						e_dep_temp.append(e_dep_temp_old)
						event_id_tot.append(event_id_tot_old)
						vol_id_tot.append(vol_id_tot_old)
						moth_id_tot.append(moth_id_tot_old)
						Strip_id_tot.append(Strip_id_tot_old)
						Si_id_tot.append(Si_id_tot_old)
						tray_id_tot.append(tray_id_tot_old)
						plane_id_tot.append(plane_id_tot_old)
						energy_dep_tot.append(energy_dep_tot_old)					
						pair_flag_tot.append(pair_flag_tot_old)

						if len(where_other_vol) != 0:
							vol_id_temp = vol_id_temp[where_other_vol]
							moth_id_temp = moth_id_temp[where_other_vol]
							Strip_id_temp = Strip_id_temp[where_other_vol]
							Si_id_temp = Si_id_temp[where_other_vol]
							tray_id_temp = tray_id_temp[where_other_vol]
							plane_id_temp = plane_id_temp[where_other_vol]
							energy_dep_temp = energy_dep_temp[where_other_vol]
							pair_flag_temp = pair_flag_temp[where_other_vol]
						else:							
							break


					where_vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]) & (Si_id_temp == 1))	
					where_vol_eq = where_vol_eq[0]
					
					where_other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]) | (Si_id_temp != 1))
					where_other_vol = where_other_vol[0]
					
					if len(where_vol_eq) != 0:
						e_dep_temp_old = np.sum(energy_dep_temp[where_vol_eq])
						event_id_tot_old = event_id_tot_temp[j]
						vol_id_tot_old = vol_id_temp[r]
						moth_id_tot_old = moth_id_temp[r]
						Si_id_tot_old = 1
						Strip_id_tot_old = Strip_id_temp[r]
						tray_id_tot_old = tray_id_temp[r]
						plane_id_tot_old = plane_id_temp[r]
						energy_dep_tot_old = e_dep_temp_old
						pair_flag_tot_old = pair_flag_temp[r]

						e_dep_temp.append(e_dep_temp_old)
						event_id_tot.append(event_id_tot_old)
						vol_id_tot.append(vol_id_tot_old)
						moth_id_tot.append(moth_id_tot_old)
						Strip_id_tot.append(Strip_id_tot_old)
						Si_id_tot.append(Si_id_tot_old)
						tray_id_tot.append(tray_id_tot_old)
						plane_id_tot.append(plane_id_tot_old)
						energy_dep_tot.append(energy_dep_tot_old)					
						pair_flag_tot.append(pair_flag_tot_old)

						if len(where_other_vol) != 0:
							vol_id_temp = vol_id_temp[where_other_vol]
							moth_id_temp = moth_id_temp[where_other_vol]
							Strip_id_temp = Strip_id_temp[where_other_vol]
							Si_id_temp = Si_id_temp[where_other_vol]
							tray_id_temp = tray_id_temp[where_other_vol]
							plane_id_temp = plane_id_temp[where_other_vol]
							energy_dep_temp = energy_dep_temp[where_other_vol]
							pair_flag_temp = pair_flag_temp[where_other_vol]
						else:
							break
			
				N_event_eq = len(where_event_eq)                                          
				if where_event_eq[N_event_eq-1] < len(event_id_tot_temp)-1:
					j = where_event_eq[N_event_eq-1]+1					
				else:
					break
				


			# apply the energy thresold
			e_dep_temp = np.array(e_dep_temp) 		
			event_id_tot = np.array(event_id_tot)
			vol_id_tot = np.array(vol_id_tot)
			moth_id_tot = np.array(moth_id_tot)
			Strip_id_tot = np.array(Strip_id_tot)
			Si_id_tot = np.array(Si_id_tot)
			tray_id_tot = np.array(tray_id_tot)
			plane_id_tot = np.array(plane_id_tot)
			energy_dep_tot = np.array(energy_dep_tot)
			pair_flag_tot = np.array(pair_flag_tot)

			
			where_eth = np.where(energy_dep_tot >= E_th)
			where_eth = where_eth[0]
			

			event_id_tot = event_id_tot[where_eth]
			vol_id_tot = vol_id_tot[where_eth]
			moth_id_tot = moth_id_tot[where_eth]
			Strip_id_tot = Strip_id_tot[where_eth]
			Si_id_tot = Si_id_tot[where_eth]
			tray_id_tot = tray_id_tot[where_eth]
			plane_id_tot = plane_id_tot[where_eth]
			energy_dep_tot = energy_dep_tot[where_eth]
			pair_flag_tot = pair_flag_tot[where_eth]


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
			

			Glob_event_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
			Glob_vol_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
			Glob_moth_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
			Glob_Strip_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
			Glob_Si_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
			Glob_tray_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
			Glob_plane_id_x_top = np.zeros((Total_vol_x_top, N_trig), dtype = np.int64)
			Glob_xpos_x_top = np.zeros((Total_vol_x_top, N_trig))
			Glob_zpos_x_top = np.zeros((Total_vol_x_top, N_trig))
			Glob_energy_dep_x_top = np.zeros((Total_vol_x_top, N_trig))
			Glob_pair_flag_x_top = np.zeros((Total_vol_x_top, N_trig))

			Glob_event_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
			Glob_vol_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
			Glob_moth_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
			Glob_Strip_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
			Glob_Si_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
			Glob_tray_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
			Glob_plane_id_y_top = np.zeros((Total_vol_y_top, N_trig), dtype = np.int64)
			Glob_ypos_y_top = np.zeros((Total_vol_y_top, N_trig))
			Glob_zpos_y_top = np.zeros((Total_vol_y_top, N_trig))
			Glob_energy_dep_y_top = np.zeros((Total_vol_y_top, N_trig))
			Glob_pair_flag_y_top = np.zeros((Total_vol_y_top, N_trig))



			for i in range(N_trig):
				Glob_vol_id_x_top[:,i] = Arch_vol_id_x_top
				Glob_moth_id_x_top[:,i] = Arch_moth_id_x_top
				Glob_Strip_id_x_top[:,i] = Arch_Strip_id_x_top
				Glob_Si_id_x_top[:,i] = Arch_Si_id_x_top
				Glob_tray_id_x_top[:,i] = Arch_tray_id_x_top
				Glob_plane_id_x_top[:,i] = Arch_plane_id_x_top
				Glob_xpos_x_top[:,i] = Arch_xpos_x_top
				Glob_zpos_x_top[:,i] = Arch_zpos_x_top
				Glob_energy_dep_x_top[:,i] = Arch_energy_dep_x_top
				Glob_pair_flag_x_top[:,i] = Arch_pair_flag_x_top

				Glob_vol_id_y_top[:,i] = Arch_vol_id_y_top
				Glob_moth_id_y_top[:,i] = Arch_moth_id_y_top
				Glob_Strip_id_y_top[:,i] = Arch_Strip_id_y_top
				Glob_Si_id_y_top[:,i] = Arch_Si_id_y_top
				Glob_tray_id_y_top[:,i] = Arch_tray_id_y_top
				Glob_plane_id_y_top[:,i] = Arch_plane_id_y_top
				Glob_ypos_y_top[:,i] = Arch_ypos_y_top
				Glob_zpos_y_top[:,i] = Arch_zpos_y_top
				Glob_energy_dep_y_top[:,i] = Arch_energy_dep_y_top
				Glob_pair_flag_y_top[:,i] = Arch_pair_flag_y_top


			
			

			j = 0
			N_ev = 0
			while 1:
				where_event_eq = np.where(event_id_tot == event_id_tot[j])
				where_event_eq = where_event_eq[0]

				
				event_id_temp = event_id_tot[where_event_eq]
				vol_id_temp = vol_id_tot[where_event_eq]
				moth_id_temp = moth_id_tot[where_event_eq]
				Strip_id_temp = Strip_id_tot[where_event_eq]
				Si_id_temp = Si_id_tot[where_event_eq]
				tray_id_temp = tray_id_tot[where_event_eq]
				plane_id_temp = plane_id_tot[where_event_eq]
				energy_dep_temp = energy_dep_tot[where_event_eq]
				pair_flag_temp = pair_flag_tot[where_event_eq]			


				# Funzione di ordinamento 1
										
				vol_sort_arr = np.argsort(vol_id_temp, kind='mergesort')
								
				vol_id_temp = vol_id_temp[vol_sort_arr]
				moth_id_temp = moth_id_temp[vol_sort_arr]
				Strip_id_temp = Strip_id_temp[vol_sort_arr]
				Si_id_temp = Si_id_temp[vol_sort_arr]
				tray_id_temp = tray_id_temp[vol_sort_arr]
				plane_id_temp = plane_id_temp[vol_sort_arr]
				energy_dep_temp = energy_dep_temp[vol_sort_arr]
				pair_flag_temp = pair_flag_temp[vol_sort_arr]				
				
				
				for z in range(Total_vol_x_top):
					where_hit_x_top_arr = np.where((Si_id_temp == 0) & (vol_id_temp == Glob_vol_id_x_top[z, N_ev]) & (moth_id_temp == Glob_moth_id_x_top[z, N_ev]))
					where_hit_x_top = where_hit_x_top_arr[0]

						
					if len(where_hit_x_top) != 0: 	
						Glob_energy_dep_x_top[z, N_ev] = energy_dep_temp[where_hit_x_top]
						Glob_pair_flag_x_top[z, N_ev] = pair_flag_temp[where_hit_x_top]
						
					where_hit_y_top_old = np.where((Si_id_temp == 1) & (vol_id_temp == Glob_vol_id_y_top[z, N_ev]) & (moth_id_temp == Glob_moth_id_y_top[z, N_ev]))
					where_hit_y_top = where_hit_y_top_old[0]

				

					if len(where_hit_y_top) != 0:
						Glob_energy_dep_y_top[z, N_ev] = energy_dep_temp[where_hit_y_top]
						Glob_pair_flag_y_top[z, N_ev] = pair_flag_temp[where_hit_y_top]

					

				N_event_eq = len(where_event_eq)
				if where_event_eq[N_event_eq-1] < len(event_id_tot)-1:
					j = where_event_eq[N_event_eq-1]+1
					N_ev = N_ev + 1
				else:
					break
			

						


			
		
			print('N_ev: '+ str(N_ev))				


			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('              Build the LEVEL 0 output            ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')		



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

			# Level 0 = energy summed
			# Level 0 = the events are sorted in tray, and Y before X within the same tray
			# energy threshold applied

			col1 = fits.Column(name='EVT_ID', format='I', array=Glob_event_id_test)	
			col2 = fits.Column(name='VOL_ID', format='J', array=Glob_vol_id_test)
			col3 = fits.Column(name='MOTH_ID', format='J', array=Glob_moth_id_test)
			col4 = fits.Column(name='TRAY_ID', format='I', array=Glob_tray_id_test)
			col5 = fits.Column(name='PLANE_ID', format='I', array=Glob_plane_id_test)
			col6 = fits.Column(name='TRK_FLAG', format='I', array=Glob_Si_id_test)
			col7 = fits.Column(name='STRIP_ID', format='J', array=Glob_Strip_id_test)
			col8 = fits.Column(name='POS', format='F20.5', array=Glob_pos_test)
			col9 = fits.Column(name='ZPOS', format='F20.5', array=Glob_zpos_test)
			col10 = fits.Column(name='E_DEP', format='F20.5', array=Glob_energy_dep_test)
			col11 = fits.Column(name='PAIR_FLAG', format='I', array=Glob_pair_flag_test)
		
			cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11])
			tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
			if os.path.exists(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+str(ifile)+'.fits'):
				os.remove(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+str(ifile)+'.fits')
				tbhdu.writeto(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
			else:
				tbhdu.writeto(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')

			fits.setval(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Creator          = Giovanni Giannella & Simone Guidotti', ext=1)

			fits.setval(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='THELSIM release  = eASTROGAM '+astrogam_version, ext=1)
		
			fits.setval(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)
		
			fits.setval(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='N_trig ='+str(N_trig), ext=1)

			fits.setval(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='ENERGY ='+ene_type, ext=1)

			fits.setval(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)
		
			fits.setval(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)
		
			fits.setval(outdir+'/L0.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = keV', ext=1)



			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - strip')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

			if os.path.exists(outdir+'/'+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
				os.remove(outdir+'/'+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
				data = open(outdir+'/'+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
			else:
				data = open(outdir+'/'+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')


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

			
			j=0
			while 1:
				where_event_eq = np.where(Glob_event_id_test == Glob_event_id_test[j])
				where_event_eq = where_event_eq[0]
				
				Glob_Si_id_test_temp = Glob_Si_id_test[where_event_eq]
				Glob_tray_id_test_temp  = Glob_tray_id_test[where_event_eq]
				Glob_plane_id_test_temp  = Glob_plane_id_test[where_event_eq]
				Glob_Strip_id_test_temp = Glob_Strip_id_test[where_event_eq]
				Glob_pos_test_temp = Glob_pos_test[where_event_eq]
				Glob_zpos_test_temp = Glob_zpos_test[where_event_eq]
				Glob_energy_dep_test_temp = Glob_energy_dep_test[where_event_eq]
				Glob_pair_flag_test_temp = Glob_pair_flag_test[where_event_eq]

				# ------------------------------------
				
				# X VIEW

				r = 0

				where_x = np.where(Glob_Si_id_test_temp == 0)
				where_x = where_x[0]				
				
				if len(where_x) != 0:				
					while r < len(where_x):
						data.write('{:d}\t'.format(Glob_event_id_test[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(ene_type))
						data.write('{:d}\t'.format(Glob_plane_id_test_temp[where_x[r]]))
						data.write('{:f}\t'.format(Glob_zpos_test_temp[where_x[r]]))
						data.write('{:d}\t'.format(0))
						data.write('{:d}\t'.format(Glob_Strip_id_test_temp[where_x[r]]))
						data.write('{:f}\t'.format(Glob_pos_test_temp[where_x[r]]))
						data.write('{:f}\t'.format(Glob_energy_dep_test_temp[where_x[r]]))
						data.write('{:d}\n'.format(Glob_pair_flag_test_temp[where_x[r]]))

						r = r + 1
				# ------------------------------------

				# Y VIEW

				r = 0

				where_y = np.where(Glob_Si_id_test_temp == 1)
				where_y = where_y[0]				
				
				if len(where_y) != 0:				
					while r < len(where_y):
						data.write('{:d}\t'.format(Glob_event_id_test[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(ene_type))
						data.write('{:d}\t'.format(Glob_plane_id_test_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_zpos_test_temp[where_y[r]]))
						data.write('{:d}\t'.format(1))
						data.write('{:d}\t'.format(Glob_Strip_id_test_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_pos_test_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_energy_dep_test_temp[where_y[r]]))
						data.write('{:d}\n'.format(Glob_pair_flag_test_temp[where_y[r]]))

						r = r + 1


				N_event_eq = len(where_event_eq)
				if where_event_eq[N_event_eq-1] < len(Glob_event_id_test)-1:
					j = where_event_eq[N_event_eq-1]+1						
				else:
					break 

			data.close()

			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('       L0.5 - cluster baricenter ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


			print('N_trig: '+ str(N_trig))


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


			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('             L0 - X-Y layers merging ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


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


			# Level 0.5 = energy summed, MIP threshold applied, strip position used

			col1 = fits.Column(name='EVT_ID', format='J', array=Glob_event_id_cluster)	
			col2 = fits.Column(name='TRAY_ID', format='I', array=Glob_tray_id_cluster)
			col3 = fits.Column(name='PLANE_ID', format='I', array=Glob_plane_id_cluster)
			col4 = fits.Column(name='TRK_FLAG', format='I', array=Glob_Si_id_cluster)
			col5 = fits.Column(name='POS', format='F.20.5', array=Glob_pos_cluster)
			col6 = fits.Column(name='ZPOS', format='F20.5', array=Glob_zpos_cluster)
			col7 = fits.Column(name='E_DEP', format='F20.5', array=Glob_energy_dep_cluster)
			col8 = fits.Column(name='PAIR_FLAG', format='I', array=Glob_pair_flag_cluster)
							
			cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8])
			tbhdu = fits.BinTableHDU.from_columns(cols)			
			
			
			if os.path.exists(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):
				os.remove(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str	(phi_type)+'.'+pol_string+str(ifile)+'.fits')
				tbhdu.writeto(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
			else:
				tbhdu.writeto(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
			
			fits.setval(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Creator          = Giovanni Giannella & Simone Guidotti', ext=1)
			
			fits.setval(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='THELSIM release  = eASTROGAM '+astrogam_version, ext=1)
			
			fits.setval(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in ='+str(N_in)+'   /Number of simulated particles', ext=1)
			
			fits.setval(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='N_trig ='+str(N_trig)+'   /Number of triggering events', ext=1)
			
			fits.setval(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='ENERGY ='+ene_type+'   /Simulated input energy', ext=1)
			
			fits.setval(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta ='+str(theta_type)+'   /Simulated input theta angle', ext=1)
			
			fits.setval(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi ='+str(phi_type)+'   /Simulated input phi angle', ext=1)
			
			fits.setval(outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = keV', ext=1)



			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - cluster')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			
			if os.path.exists(outdir+'/'+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
				os.remove(outdir+'/'+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
				data = open(outdir+'/'+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
			else:
				data = open(outdir+'/'+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')


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

			totalstrips_before = 0
			j=0
			while 1:
				where_event_eq = np.where(Glob_event_id_cluster == Glob_event_id_cluster[j])
				where_event_eq = where_event_eq[0] 

				Glob_Si_id_cluster_temp = Glob_Si_id_cluster[where_event_eq]
				Glob_tray_id_cluster_temp  = Glob_tray_id_cluster[where_event_eq]
				Glob_plane_id_cluster_temp  = Glob_plane_id_cluster[where_event_eq]
				Glob_energy_dep_cluster_temp = Glob_energy_dep_cluster[where_event_eq]
				Glob_pos_cluster_temp = Glob_pos_cluster[where_event_eq]
				Glob_zpos_cluster_temp = Glob_zpos_cluster[where_event_eq]
				Glob_pair_flag_cluster_temp = Glob_pair_flag_cluster[where_event_eq]
				Glob_Strip_number_cluster_temp = Glob_Strip_number_cluster[where_event_eq]


				# ------------------------------------
				
				# X VIEW

				r = 0

				where_x = np.where(Glob_Si_id_cluster_temp == 0)
				where_x = where_x[0]				
				
				if len(where_x) != 0:				
					while r < len(where_x):
						data.write('{:d}\t'.format(Glob_event_id_cluster[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(ene_type))
						data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_x[r]]))
						data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_x[r]]))
						data.write('{:d}\t'.format(0))
						data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_x[r]]))
						data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_x[r]]))
						data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_x[r]]))
						data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_x[r]]))

						r = r + 1
				# ------------------------------------

				# Y VIEW

				r = 0

				where_y = np.where(Glob_Si_id_cluster_temp == 1)
				where_y = where_y[0]				
				
				if len(where_y) != 0:				
					while r < len(where_y):
						data.write('{:d}\t'.format(Glob_event_id_cluster[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(ene_type))
						data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_y[r]]))
						data.write('{:d}\t'.format(1))
						data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_y[r]]))
						data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_y[r]]))
						data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_y[r]]))

						r = r + 1


				N_event_eq = len(where_event_eq)
				if where_event_eq[N_event_eq-1] < len(Glob_event_id_cluster)-1:
					j = where_event_eq[N_event_eq-1]+1						
				else:
					break 


			data.close()
					
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - pairs')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			
			if os.path.exists(outdir+'/'+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
				os.remove(outdir+'/'+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
				data = open(outdir+'/'+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
			else:
				data = open(outdir+'/'+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')


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

			totalstrips_before = 0
			j=0
			while 1:
				where_event_eq = np.where(Glob_event_id_cluster == Glob_event_id_cluster[j])
				where_event_eq = where_event_eq[0] 

				Glob_Si_id_cluster_temp = Glob_Si_id_cluster[where_event_eq]
				Glob_tray_id_cluster_temp  = Glob_tray_id_cluster[where_event_eq]
				Glob_plane_id_cluster_temp  = Glob_plane_id_cluster[where_event_eq]
				Glob_energy_dep_cluster_temp = Glob_energy_dep_cluster[where_event_eq]
				Glob_pos_cluster_temp = Glob_pos_cluster[where_event_eq]
				Glob_zpos_cluster_temp = Glob_zpos_cluster[where_event_eq]
				Glob_pair_flag_cluster_temp = Glob_pair_flag_cluster[where_event_eq]
				Glob_Strip_number_cluster_temp = Glob_Strip_number_cluster[where_event_eq]


				# ------------------------------------
				
				# X VIEW

				r = 0

				where_x = np.where(Glob_Si_id_cluster_temp == 0)
				where_x = where_x[0]				
				
				if len(where_x) != 0:				
					while r < len(where_x):
						if Glob_pair_flag_cluster_temp[where_x[r]] == 1:
							data.write('{:d}\t'.format(Glob_event_id_cluster[j]))
							data.write('{:d}\t'.format(theta_type))
							data.write('{:d}\t'.format(phi_type))
							data.write('{:s}\t'.format(ene_type))
							data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_x[r]]))
							data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_x[r]]))
							data.write('{:d}\t'.format(0))
							data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_x[r]]))
							data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_x[r]]))
							data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_x[r]]))
							data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_x[r]]))

						r = r + 1
				# ------------------------------------

				# Y VIEW

				r = 0

				where_y = np.where(Glob_Si_id_cluster_temp == 1)
				where_y = where_y[0]				
				
				if len(where_y) != 0:				
					while r < len(where_y):
						if Glob_pair_flag_cluster_temp[where_y[r]] == 1:
							data.write('{:d}\t'.format(Glob_event_id_cluster[j]))
							data.write('{:d}\t'.format(theta_type))
							data.write('{:d}\t'.format(phi_type))
							data.write('{:s}\t'.format(ene_type))
							data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_y[r]]))
							data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_y[r]]))
							data.write('{:d}\t'.format(1))
							data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_y[r]]))
							data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_y[r]]))
							data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_y[r]]))
							data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_y[r]]))

						r = r + 1


				N_event_eq = len(where_event_eq)
				if where_event_eq[N_event_eq-1] < len(Glob_event_id_cluster)-1:
					j = where_event_eq[N_event_eq-1]+1						
				else:
					break 					

			data.close()

			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - compton')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
								   
			if os.path.exists(outdir+'/'+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
				os.remove(outdir+'/'+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
				data = open(outdir+'/'+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
			else:
				data = open(outdir+'/'+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
									   
									   
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
									   
			totalstrips_before = 0
			j=0
			while 1:
				where_event_eq = np.where(Glob_event_id_cluster == Glob_event_id_cluster[j])
				where_event_eq = where_event_eq[0] 
									   
				Glob_Si_id_cluster_temp = Glob_Si_id_cluster[where_event_eq]
				Glob_tray_id_cluster_temp  = Glob_tray_id_cluster[where_event_eq]
				Glob_plane_id_cluster_temp  = Glob_plane_id_cluster[where_event_eq]
				Glob_energy_dep_cluster_temp = Glob_energy_dep_cluster[where_event_eq]
				Glob_pos_cluster_temp = Glob_pos_cluster[where_event_eq]
				Glob_zpos_cluster_temp = Glob_zpos_cluster[where_event_eq]
				Glob_pair_flag_cluster_temp = Glob_pair_flag_cluster[where_event_eq]
				Glob_Strip_number_cluster_temp = Glob_Strip_number_cluster[where_event_eq]
									   
									   
				# ------------------------------------
									   
				# X VIEW
									   
				r = 0

				where_x = np.where(Glob_Si_id_cluster_temp == 0)
				where_x = where_x[0]				
									   
				if len(where_x) != 0:
					while r < len(where_x):
						if Glob_pair_flag_cluster_temp[where_x[r]] == 2:
							data.write('{:d}\t'.format(Glob_event_id_cluster[j]))
							data.write('{:d}\t'.format(theta_type))
							data.write('{:d}\t'.format(phi_type))
							data.write('{:s}\t'.format(ene_type))
							data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_x[r]]))
							data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_x[r]]))
							data.write('{:d}\t'.format(0))
							data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_x[r]]))
							data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_x[r]]))
							data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_x[r]]))
							data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_x[r]]))
											  
						r = r + 1
							
				# ------------------------------------
											
				# Y VIEW
											  
				r = 0
											  
				where_y = np.where(Glob_Si_id_cluster_temp == 1)
				where_y = where_y[0]				
											  
				if len(where_y) != 0:				
					while r < len(where_y):
						if Glob_pair_flag_cluster_temp[where_y[r]] == 2:
							data.write('{:d}\t'.format(Glob_event_id_cluster[j]))
							data.write('{:d}\t'.format(theta_type))
							data.write('{:d}\t'.format(phi_type))
							data.write('{:s}\t'.format(ene_type))
							data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_y[r]]))
							data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_y[r]]))
							data.write('{:d}\t'.format(1))
							data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_y[r]]))
							data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_y[r]]))
							data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_y[r]]))
							data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_y[r]]))
														 
						r = r + 1
														 
														 
				N_event_eq = len(where_event_eq)
				if where_event_eq[N_event_eq-1] < len(Glob_event_id_cluster)-1:
					j = where_event_eq[N_event_eq-1]+1						
				else:
					break 					


			data.close()

	if cal_flag == 1:
		
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
			#print(where_event_eq[N_event_eq-1])
			if where_event_eq[N_event_eq-1] < len(event_id_cal)-1:
				j = where_event_eq[N_event_eq-1]+1
			else:
				break
		
		
		event_id_tot_cal = np.array((event_id_tot_cal), dtype = np.int64)
		bar_id_tot = np.array((bar_id_tot), dtype = np.int64)
		bar_ene_tot = np.array(bar_ene_tot)
		
		col1 = fits.Column(name='EVT_ID', format='I', array=event_id_tot_cal)	
		col2 = fits.Column(name='BAR_ID', format='I', array=bar_id_tot)
		col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=bar_ene_tot)

		
		cols = fits.ColDefs([col1,col2,col3])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):
			os.remove(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
			tbhdu.writeto(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')

		fits.setval(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

		fits.setval(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)



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
		
		col1 = fits.Column(name='EVT_ID', format='I', array=event_id_tot_cal_sum)	
		col2 = fits.Column(name='BAR_ENERGY', format='F20.15', array=bar_ene_tot_sum)

		cols = fits.ColDefs([col1,col2])
		tbhdu = fits.BinTableHDU.from_columns(cols)


		if os.path.exists(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):
			os.remove(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
			tbhdu.writeto(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')

		fits.setval(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

		fits.setval(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)

		

	if ac_flag == 1:

		event_id_ac = np.array(event_id_ac)
		vol_id_ac = np.array(vol_id_ac)
		en_dep_ac = np.array(en_dep_ac)
		moth_id_ac =np.array(moth_id_ac)


		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('                          AC')
		print('                  Summing the energy                ')
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

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
					AC_panel_old[j] = 'S'
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


		col1 = fits.Column(name='EVT_ID', format='I', array=event_id_tot_ac)	
		col2 = fits.Column(name='AC_PANEL', format='A', array=AC_panel)
		col3 = fits.Column(name='AC_SUBPANEL', format='I', array=AC_subpanel)
		col4 = fits.Column(name='E_DEP', format='F20.15', array=energy_dep_tot_ac)
		
		cols = fits.ColDefs([col1,col2,col3,col4])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits'):
			os.remove(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
			tbhdu.writeto(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits')


		fits.setval(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+ene_type, ext=1)

		fits.setval(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(outdir+'/G4.AC.eASTROGAMM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)


############    passo a file fits successivo   ####################

	t.close()

	end = time.time() - start

	print('Time: ' + str(end))
				
	ifile = ifile + 1



	













