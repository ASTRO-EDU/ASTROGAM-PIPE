from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import os

##############################
#     parametri iniziali     #
##############################

astrogam_version = 'V1.0'   # Enter eASTROGAM release (e.g. V1.0):
bogemms_tag = 211         # Enter BoGEMMS release (e.g. 211):
sim_type = 0            # Enter simulation type [0 = Mono, 1 = Range, 2 = Chen, 3: Vela, 4: Crab, 4: G400]:
py_list = 400          # Enter the Physics List [0 = QGSP_BERT_EMV, 100 = ARGO, 300 = FERMI, 400 = ASTROMEV]:
N_in = 100000                # Enter the number of emitted particles:
part_type = "ph"           # Enter the particle type [ph = photons, mu = muons, g = geantino, p = proton, el = electron]:
n_fits = 2              # Enter number of FITS files:
ene_range = 3          # Enter energy distribution [0 = MONO, 1 = POW, 2 = EXP, 3 = LIN]:
ene_min = .5             # Enter miminum energy [MeV]:
ene_max = 100            # Enter maximum energy [MeV]:
ang_type = "UNI"           # Enter the angular distribution [e.g. UNI, ISO]:
theta_type = 30         # Enter theta:
phi_type = 225           # Enter phi:
pol_type = 1           # Is the source polarized? [0 = false, 1 = true]: 
pol_angle = 20          # Enter Polarization angle:
source_g = 0           # Enter source geometry [0 = Point, 1 = Plane]:
isStrip = 1            # Strip/Pixels activated?:
repli = 1              # Strips/Pixels replicated?:
cal_flag = 1           # Is Cal present? [0 = false, 1 = true]:
ac_flag = 0            # Is AC present? [0 = false, 1 = true]:
passive_flag = 0       # Is Passive present? [0 = false, 1 = true]:
energy_thresh = 15      # Enter energy threshold [keV]:

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

#print(N_tray, panel_S, strip_side, energy_thresh)


#parte lettura file fits

n_fits = os.listdir("/home/gianni/eASTROGAMSimAnalysis_STUDENT")
ifile = 0
print(n_fits)
while ifile < len(n_fits):
	t = fits.open(n_fits[ifile])

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

	i = 0
	max_dim = len(tbdata)

	while i < max_dim:
		
		vol_id = volume_id[i]
        	moth_id = mother_id[i]
		energy_dep = e_dep[i]
		#  Reading the tracker (events with E > 0)
       		        
		if vol_id >= tracker_bottom_vol_start or moth_id >= tracker_bottom_vol_start:

			event_id = evt_id[i]
        		track_id = trk_id[i]
       			vol_name = volume_name[i]
			en_dep = e_dep[i]
			x_en = x_ent[i]
	        	y_en = y_ent[i]
	        	z_en = z_ent[i]
	        	x_ex = x_exit[i]
	        	y_ex = y_exit[i]
	        	z_ex = z_exit[i]
	        	e_kin_en = e_kin_ent[i]
	        	e_kin_ex = e_kin_exit[i]
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
			
			if part_type == 'g':
				e_dep = 100.
				en_dep = e_dep[i]
			elif energy_dep > 0.:
				continue

			theta_ent = (180./math.pi)*math.acos(-(cos_x_angle_ent))
          		phi_ent = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))

          		theta_exit = (180./math.pi)*math.acos(-(cos_z_angle_exit))
          		phi_exit = (180./math.pi)*math.atan((cos__angle_exit)/(cos_x_angle_exit))

          		child_id = parent_trk_id[i]
          		proc_id = process_id[i]


		# Reading the Calorimeter (controllare separazione geantino e altro)

		if cal_flag == 1:
			if vol_id >= cal_vol_start and vol_id <= cal_vol_end:
				if part_type == 'g' or energy_dep > 0.:
				event_id_cal = evt_id[i]
	        		track_id_cal = trk_id[i]
				vol_id_cal = vol_id[i]	       			
				vol_name_cal = volume_name[i]
				en_dep_cal = e_dep[i]
				x_en_cal = x_ent[i]
		        	y_en_cal = y_ent[i]
		        	z_en_cal = z_ent[i]
		        	x_ex_cal = x_exit[i]
		        	y_ex_cal = y_exit[i]
		        	z_ex_cal = z_exit[i]
		        	e_kin_en_cal = e_kin_ent[i]
		        	e_kin_ex_cal = e_kin_exit[i]
		        	cos_x_angle_ent_cal = mdx_ent[i]
		        	cos_y_angle_ent_cal = mdy_ent[i]
		        	cos_z_angle_ent_cal = mdz_ent[i]
		        	cos_x_angle_exit_cal = mdx_exit[i]
		        	cos_y_angle_exit_cal = mdy_exit[i]
		        	cos_z_angle_exit_cal = mdz_exit[i]
		        	gtime_en_cal = gtime_ent[i]
		        	gtime_ex_cal = gtime_exit[i]
		        	part_id_cal = particle_id[i]
		        	part_name_cal = particle_name[i]
		        	proc_name_cal = process_name[i]
				en_dep_cal = e_dep[i]
				
				theta_ent_cal = (180./math.pi)*math.acos(-(cos_x_angle_ent))
	          		phi_ent_cal = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
	
	          		theta_exit_cal = (180./math.pi)*math.acos(-(cos_z_angle_exit))
	          		phi_exit_cal = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

	          		child_id_cal = parent_trk_id[i]
	          		proc_id_cal = process_id[i]


		# Reading the AC

		if ac_flag == 1:
			if vol_id >= ac_vol_start and vol_id <= ac_vol_end:
				if part_type == 'g' or energy_dep > 0.:
					event_id_ac = evt_id[i]
					vol_id_cal = vol_id[i]

					energy_dep_ac = e_dep[i]

					ent_x_ac = x_ent[i]
        				ent_y_ac = y_ent[i]
        				ent_z_ac = z_ent[i]
        				exit_x_ac = x_exit[i]
        				exit_y_ac = y_exit[i]
					exit_z_ac = z_exit[i]

					theta_ent_ac = (180./math.pi)*math.acos(-(cos_x_angle_ent))
	          			phi_ent_ac = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
	
	          			theta_exit_ac = (180./math.pi)*math.acos(-(cos_z_angle_exit))
	          			phi_exit_ac = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

	          			child_id_ac = parent_trk_id[i]
	          			proc_id_ac = process_id[i]

						if isStrip == 1:
							moth_id_ac = mother_id[i]
						else:
							moth_id_ac = 0







		i = i + 1	









	ifile = ifile + 1















