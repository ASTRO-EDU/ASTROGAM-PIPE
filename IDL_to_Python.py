from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import pickle
from astropy.table import Table

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

n_fits = os.listdir('./eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+dir_cal+dir_passive+'/'+ene_type+'MeV.'+sim_name+'.'+str(theta_type)+'theta.'+pol_string+str(N_in)+part_type)

filepath = '/eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+dir_cal+dir_passive+'/'+ene_type+'MeV.'+sim_name+'.'+str(theta_type)+'theta.'+pol_string+str(N_in)+part_type

print('eASTROGAM simulation path:' + filepath)

outdir = ('./eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+str(sim_name)+'/'+str(ene_type)+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV')

if not os.path.exists('./eASTROGAM'+astrogam_version+sdir+'/theta'+str(theta_type)+'/'+stripDir+py_dir+'/'+str(sim_name)+'/'+str(ene_type)+'MeV/'+str(N_in)+part_type+dir_cal+dir_passive+'/'+str(energy_thresh)+'keV'):
	out_dir = os.makedirs(outdir,0777)

ifile = 0
while ifile < len(n_fits):

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

	i = 0
	max_dim = len(tbdata)
	
	vol_id_track = []	
	moth_id_track = []
	event_id_track = []
	en_dep_track = []
	x_en_track = []
	y_en_track = []
	z_en_track = []
	x_ex_track = []
	y_ex_track = []
	z_ex_track = []
	child_id_track = []
	proc_id_track = []
	theta_ent_track = []
	phi_ent_track = []
	theta_exit_track = []
	phi_exit_track = []	
	x_tr = []
	y_tr = []
	z_tr = []



	while i < max_dim:
		
		vol_id = volume_id[i]         #volume per condizioni
		moth_id = mother_id[i]
	        energy_dep = e_dep[i]         #energia per condizioni
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
				event_id = evt_id[i]
				track_id = trk_id[i]
				vol_name = volume_name[i]
				en_dep = e_dep[i]
				x_en = (x_ent[i])/10.
				y_en = (y_ent[i])/10.
				z_en = (z_ent[i])/10.
				x_ex = (x_exit[i])/10.
				y_ex = (y_exit[i])/10.
				z_ex = (z_exit[i])/10.
				e_kin_en = e_kin_ent[i]
				e_kin_ex = e_kin_exit[i]
				v_id_track = volume_id[i]
				m_id_track = mother_id[i]
				e_dep = 100.
				en_dep = e_dep[i]    #controllare

				theta_ent = (180./math.pi)*math.acos(-(cos_x_angle_ent))
				phi_ent = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
	
				theta_exit = (180./math.pi)*math.acos(-(cos_z_angle_exit))
				phi_exit = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

				child_id = parent_trk_id[i]
				proc_id = process_id[i]
			
				x_position = x_en + ((x_ex - x_en)/2.)
				y_position = y_en + ((y_ex - y_en)/2.)
				z_position = z_en + ((z_ex - z_en)/2.)

				vol_id_track.append(v_id_track)
				moth_id_track.append(m_id_track)
				event_id_track.append(event_id)
				en_dep_track.append(en_dep)			
				x_en_track.append(x_en)
				y_en_track.append(y_en)
				z_en_track.append(z_en)
				x_ex_track.append(x_ex)
				y_ex_track.append(x_ex)
				z_ex_track.append(x_ex)
				child_id_track.append(child_id)
				proc_id_track.append(proc_id)
				theta_ent_track.append(theta_ent)
				phi_ent_track.append(phi_ent)
				theta_exit_track.append(theta_exit)
				phi_exit_track.append(phi_exit)
				x_tr.append(x_position)
				y_tr.append(y_position)
				z_tr.append(z_position)

			if energy_dep > 0.:

				event_id = evt_id[i]
				track_id = trk_id[i]
				vol_name = volume_name[i]
				en_dep = e_dep[i]
				x_en = (x_ent[i])/10.
				y_en = (y_ent[i])/10.
				z_en = (z_ent[i])/10.
				x_ex = (x_exit[i])/10.
				y_ex = (y_exit[i])/10.
				z_ex = (z_exit[i])/10.
				e_kin_en = e_kin_ent[i]
				e_kin_ex = e_kin_exit[i]
				v_id_track = volume_id[i]
				m_id_track = mother_id[i]
				en_dep = e_dep[i]    
				
				theta_ent = (180./math.pi)*math.acos(-(cos_x_angle_ent))
				phi_ent = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))

				theta_exit = (180./math.pi)*math.acos(-(cos_z_angle_exit))
				phi_exit = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

				child_id = parent_trk_id[i]
				proc_id = process_id[i]
			
				x_position = x_en + ((x_ex - x_en)/2.)
				y_position = y_en + ((y_ex - y_en)/2.)
				z_position = z_en + ((z_ex - z_en)/2.)

				vol_id_track.append(v_id_track)
				moth_id_track.append(m_id_track)
				event_id_track.append(event_id)
				en_dep_track.append(en_dep)			
				x_en_track.append(x_en)
				y_en_track.append(y_en)
				z_en_track.append(z_en)
				x_ex_track.append(x_ex)
				y_ex_track.append(x_ex)
				z_ex_track.append(x_ex)
				child_id_track.append(child_id)
				proc_id_track.append(proc_id)
				theta_ent_track.append(theta_ent)
				phi_ent_track.append(phi_ent)
				theta_exit_track.append(theta_exit)
				phi_exit_track.append(phi_exit)
				x_tr.append(x_position)
				y_tr.append(y_position)
				z_tr.append(z_position)
				

		# Reading the Calorimeter (controllare separazione geantino e altro)
		if cal_flag == 1:
            #print(i,'a', energy_dep)
			if vol_id >= cal_vol_start and vol_id <= cal_vol_end:
                #print('b', energy_dep)
	                	if part_type == 'g' or energy_dep > 0.:
                    #print('c')
	                    		event_id_cal = evt_id[i]
	                    		track_id_cal = trk_id[i]
	                    		vol_id_cal = volume_id[i]
	                    		vol_name_cal = volume_name[i]
		                  	en_dep_cal = e_dep[i]
	                    		x_en_cal = (x_ent[i])/10.
	                    		y_en_cal = (y_ent[i])/10.
	                    		z_en_cal = (z_ent[i])/10.
	                    		x_ex_cal = (x_exit[i])/10.
	                    		y_ex_cal = (y_exit[i])/10.
	                    		z_ex_cal = (z_exit[i])/10.
	                    		e_kin_en_cal = e_kin_ent[i]
	                    		e_kin_ex_cal = e_kin_exit[i]
	                    		en_dep_cal = e_dep[i]
				
	                    		theta_ent_cal = (180./math.pi)*math.acos(-(cos_x_angle_ent))
	                    		phi_ent_cal = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
		
	                    		theta_exit_cal = (180./math.pi)*math.acos(-(cos_z_angle_exit))
	                    		phi_exit_cal = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

	                    		child_id_cal = parent_trk_id[i]
	                    		proc_id_cal = process_id[i]
	                    		#print(phi_exit_cal)

		# Reading the AC

		if ac_flag == 1:
			if vol_id >= ac_vol_start and vol_id <= ac_vol_end:
				if part_type == 'g' or energy_dep > 0.:
					event_id_ac = evt_id[i]
					vol_id_cal = vol_id[i]

					energy_dep_ac = e_dep[i]

					ent_x_ac = (x_ent[i])/10.
					ent_y_ac = (y_ent[i])/10.
					ent_z_ac = (z_ent[i])/10.
					exit_x_ac = (x_exit[i])/10.
					exit_y_ac = (y_exit[i])/10.
					exit_z_ac = (z_exit[i])/10.

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
		max_dim = len(vol_id_track)
		while j < max_dim:
			if isStrip == 1:
				if repli == 1:
					Strip_y = vol_id_track[j]
					tray = moth_id_track[j]/tracker_bottom_vol_start
					invert_tray_id = (N_tray - tray)+1
					vol_id_temp = moth_id_track[j] - (tracker_bottom_vol_start*tray + tracker_top_bot_diff) # removing 1000000xn_tray + 90000					#
					Strip_x = vol_id_temp
					plane = invert_tray_id	

					plane_id.append(plane)
					Strip_id_y.append(Strip_y)
					Strip_id_x.append(Strip_x)
					tray_id.append(tray)

			else:	
					

				Strip_y = 0
				tray = vol_id_track[j]/tracker_bottom_vol_start
				invert_tray_id = (N_tray - tray)+1
				Strip_x= 0
				plane = invert_tray_id				
				
				plane_id.append(plane)
				Strip_id_y.append(Strip_y)
				Strip_id_x.append(Strip_x)
				tray_id.append(tray)

			j = j+1
		#print(Strip_id_x)


	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('                             Tracker   '                     )
	print('           Saving the Tracker raw hits (fits and .dat)      ')
	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	if astrogam_version == 'V1.0':
 	
		event_id_tr = np.array(event_id_track)
		vol_id_tr = np.array(vol_id_track)
		moth_id_tr = np.array(moth_id_track)
		tr_id = np.array(tray_id)
		pl_id = np.array(plane_id)
		Str_id_x = np.array(Strip_id_x)
		Str_id_y = np.array(Strip_id_y)
		en_dep_tr = np.array(en_dep_track)
		x_en_tr = np.array(x_en_track)
		y_en_tr = np.array(y_en_track)
		z_en_tr = np.array(z_en_track)
		x_ex_tr = np.array(x_ex_track)
		y_ex_tr = np.array(y_ex_track)
		z_ex_tr = np.array(z_ex_track)	
		child_id_tr = np.array(child_id_track)
		proc_id_tr = np.array(proc_id_track)
		x_pos = np.array(x_tr)
		y_pos = np.array(y_tr)
		z_pos = np.array(z_tr)

		


		col1 = fits.Column(name='EVT_ID', format='I', array=event_id_tr)	
		col2 = fits.Column(name='VOL_ID', format='I', array=vol_id_tr)
		col3 = fits.Column(name='MOTH_ID', format='J', array=moth_id_tr)
		col4 = fits.Column(name='TRAY_ID', format='I', array=tr_id)
		col5 = fits.Column(name='PLANE_ID', format='I', array=pl_id)
		col6 = fits.Column(name='STRIP_ID_X', format='I', array=Str_id_x)
		col7 = fits.Column(name='STRIP_ID_Y', format='I', array=Str_id_y)
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
			if os.path.exists(outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
				os.remove(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
				data = open(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
			else:
				data = open(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
	
			

			j = 0
			while j < len(event_id_tr):
				event_eq = np.where(event_id_tr == event_id_tr[j])
				
				where_event_eq = event_eq[0]
				plane_id_temp = pl_id[where_event_eq]
				Cluster_x_temp  = x_pos[where_event_eq]
				Cluster_y_temp  = y_pos[where_event_eq]
				Cluster_z_temp  = z_pos[where_event_eq]
				e_dep_x_temp  = (en_dep_tr[where_event_eq])/2.
				e_dep_y_temp  = (en_dep_tr[where_event_eq])/2.
				child_temp = child_id_tr[where_event_eq]
				proc_temp = proc_id_tr[where_event_eq]
				print(plane_id_temp)
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

				j_max = max(where_event_eq)
				
				j = j_max + 1
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
			while 1:           #j < len(event_id_tr):
				event_eq = np.where(event_id_tr == event_id_tr[j])
				where_event_eq = event_eq[0]

				vol_id_temp = vol_id_tr[where_event_eq]
				moth_id_temp  = moth_id_tr[where_event_eq]
				Strip_id_x_temp  = Str_id_x[where_event_eq]
				Strip_id_y_temp  = Str_id_y[where_event_eq]
				tray_id_temp  = tr_id[where_event_eq]
				plane_id_temp  = pl_id[where_event_eq]
				energy_dep_temp = en_dep_tr[where_event_eq]
				child_id_temp = child_id_tr[where_event_eq]
				proc_id_temp = proc_id_tr[where_event_eq]
				
				r = 0												
				while 1:
					#where_vol_eq e where_other_vol numero array inferiore rispetto IDL
					vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]))										
					where_vol_eq = vol_eq[0]	
					
					other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]))
					where_other_vol = other_vol[0]
					
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

					w_pair = np.where((all_child == 1) & (all_proc == 7))
					where_pair = w_pair[0]

					if where_pair != []:
						pair_flag_tot_old = 1
					else:
						pair_flag_tot_old = 0				#il terzultimo valore di IDL e' 1 quando dovrebbe essere 0 perche nel where_pair e' -1   
					pair_flag_tot.append(pair_flag_tot_old)	
 
					w_compton = np.where((all_child == 1) & (all_proc == 3))
					where_compton = w_compton[0]

					if where_compton != []:
						pair_flag_tot_old = 2
					else:
						pair_flag_tot_old = 0
				
					pair_flag_tot.append(pair_flag_tot_old)	

					if where_other_vol != []:
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
				#j_max = max(where_event_eq)
				#j = j_max + 1
			
				N_event_eq = len(where_event_eq)                                          
				if where_event_eq[N_event_eq-1] < len(event_id_tr)-1:
					j = where_event_eq[N_event_eq-1]+1					
				else:
					break

			event_id_tot_temp = np.zeros(2*len(event_id_tot))
			vol_id_tot_temp = np.zeros(2*len(event_id_tot))
			moth_id_tot_temp = np.zeros(2*len(event_id_tot))
			Strip_id_tot_temp = np.zeros(2*len(event_id_tot))
			Si_id_tot_temp = np.zeros(2*len(event_id_tot))
			tray_id_tot_temp = np.zeros(2*len(event_id_tot))
			plane_id_tot_temp = np.zeros(2*len(event_id_tot))
			energy_dep_tot_temp = np.zeros(2*len(event_id_tot))
			pair_flag_tot_temp = np.zeros(2*len(event_id_tot))

			jev = 0
			while jev < len(event_id_tot):
				ev_index = jev * 2

				event_id_tot_temp = np.delete(event_id_tot_temp, ev_index)
				event_id_tot_temp = np.insert(event_id_tot_temp, ev_index, event_id_tot[jev])

				event_id_tot_temp = np.delete(event_id_tot_temp, ev_index + 1)
				event_id_tot_temp = np.insert(event_id_tot_temp, ev_index + 1, event_id_tot[jev])

				vol_id_tot_temp = np.delete(vol_id_tot_temp, ev_index)
				vol_id_tot_temp = np.insert(vol_id_tot_temp, ev_index, Strip_id_x_tot[jev])

				vol_id_tot_temp = np.delete(vol_id_tot_temp, ev_index + 1)
				vol_id_tot_temp = np.insert(vol_id_tot_temp, ev_index + 1, Strip_id_y_tot[jev])

				moth_id_tot_temp = np.delete(moth_id_tot_temp, ev_index)
				moth_id_tot_temp = np.insert(moth_id_tot_temp, ev_index, moth_id_tot[jev] - Strip_id_x_tot[jev])

				moth_id_tot_temp = np.delete(moth_id_tot_temp, ev_index + 1)
				moth_id_tot_temp = np.insert(moth_id_tot_temp, ev_index + 1, moth_id_tot[jev] - Strip_id_x_tot[jev])

				Strip_id_tot_temp = np.delete(Strip_id_tot_temp, ev_index)
				Strip_id_tot_temp = np.insert(Strip_id_tot_temp, ev_index, Strip_id_x_tot[jev])

				Strip_id_tot_temp = np.delete(Strip_id_tot_temp, ev_index + 1)
				Strip_id_tot_temp = np.insert(Strip_id_tot_temp, ev_index + 1, Strip_id_y_tot[jev])

				Si_id_tot_temp = np.delete(Si_id_tot_temp, ev_index)
				Si_id_tot_temp = np.insert(Si_id_tot_temp, ev_index, 0)
				
				Si_id_tot_temp = np.delete(Si_id_tot_temp, ev_index + 1)
				Si_id_tot_temp = np.insert(Si_id_tot_temp, ev_index + 1, 1)

				tray_id_tot_temp = np.delete(tray_id_tot_temp, ev_index)
				tray_id_tot_temp = np.insert(tray_id_tot_temp, ev_index, tray_id_tot[jev])

				tray_id_tot_temp = np.delete(tray_id_tot_temp, ev_index + 1)
				tray_id_tot_temp = np.insert(tray_id_tot_temp, ev_index + 1, tray_id_tot[jev])

				plane_id_tot_temp = np.delete(plane_id_tot_temp, ev_index)
				plane_id_tot_temp = np.insert(plane_id_tot_temp, ev_index, plane_id_tot[jev])

				plane_id_tot_temp = np.delete(plane_id_tot_temp, ev_index + 1)
				plane_id_tot_temp = np.insert(plane_id_tot_temp, ev_index + 1, plane_id_tot[jev])

				tray_id_tot_temp = np.delete(tray_id_tot_temp, ev_index)
				tray_id_tot_temp = np.insert(tray_id_tot_temp, ev_index, tray_id_tot[jev])

				energy_dep_tot_temp = np.delete(energy_dep_tot_temp, ev_index)
				energy_dep_tot_temp = np.insert(energy_dep_tot_temp, ev_index, energy_dep_tot[jev]/2.)

				energy_dep_tot_temp = np.delete(energy_dep_tot_temp, ev_index + 1)
				energy_dep_tot_temp = np.insert(energy_dep_tot_temp, ev_index + 1, energy_dep_tot[jev]/2.)

				pair_flag_tot_temp = np.delete(pair_flag_tot_temp, ev_index)
				pair_flag_tot_temp = np.insert(pair_flag_tot_temp, ev_index, pair_flag_tot[jev])

				pair_flag_tot_temp = np.delete(pair_flag_tot_temp, ev_index + 1)
				pair_flag_tot_temp = np.insert(pair_flag_tot_temp, ev_index + 1, pair_flag_tot[jev])

				
				jev = jev + 1
						
			#
			# Summing the energy along the strip
			#

			e_dep_temp_Si_0 = [] 		
			event_id_tot_Si_0 = []
			vol_id_tot_Si_0 = []
			moth_id_tot_Si_0 = []
			Strip_id_x_tot_Si_0 = []
			Strip_id_y_tot_Si_0 = []
			tray_id_tot_Si_0 = []
			plane_id_tot_Si_0 = []
			energy_dep_tot_Si_0 = []
			pair_flag_tot_Si_0 = []

			e_dep_temp_Si_1 = [] 		
			event_id_tot_Si_1 = []
			vol_id_tot_Si_1 = []
			moth_id_tot_Si_1 = []
			Strip_id_x_tot_Si_1 = []
			Strip_id_y_tot_Si_1 = []
			tray_id_tot_Si_1 = []
			plane_id_tot_Si_1 = []
			energy_dep_tot_Si_1 = []
			pair_flag_tot_Si_1 = []

			j = 0
			while 1:           
				event_eq = np.where(event_id_tot_temp == event_id_tot_temp[j])
				where_event_eq = event_eq[0]

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
					#where_vol_eq e where_other_vol numero array inferiore rispetto IDL
					vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]) & (Si_id_temp == 0))	
					where_vol_eq = vol_eq[0]
					
					other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]) | (Si_id_temp != 0))
					where_other_vol = other_vol[0]

					if where_vol_eq != []:
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

						e_dep_temp_Si_0.append(e_dep_temp_old)
						event_id_tot_Si_0.append(event_id_tot_old)
						vol_id_tot_Si_0.append(vol_id_tot_old)
						moth_id_tot_Si_0.append(moth_id_tot_old)
						Strip_id_x_tot_Si_0.append(Strip_id_x_tot_old)
						Strip_id_y_tot_Si_0.append(Strip_id_y_tot_old)
						tray_id_tot_Si_0.append(tray_id_tot_old)
						plane_id_tot_Si_0.append(plane_id_tot_old)
						energy_dep_tot_Si_0.append(energy_dep_tot_old)					
					

						if where_other_vol != []:
							vol_id_temp = vol_id_temp[where_other_vol]
							moth_id_temp = moth_id_temp[where_other_vol]
							Strip_id_temp = Strip_id_temp[where_other_vol]
							Si_id_temp = Si_id_temp[where_other_vol]
							tray_id_temp = tray_id_temp[where_other_vol]
							plane_id_temp = plane_id_temp[where_other_vol]
							energy_dep_temp = energy_dep_temp[where_other_vol]
							pair_flag_temp = pair_flag_temp[where_other_vol]
						else:
							#print('exit'
							break


					vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]) & (Si_id_temp == 1))	
					where_vol_eq = vol_eq[0]
					
					other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]) | (Si_id_temp != 1))
					where_other_vol = other_vol[0]
					
					if where_vol_eq != []:
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

						e_dep_temp_Si_1.append(e_dep_temp_old)
						event_id_tot_Si_1.append(event_id_tot_old)
						vol_id_tot_Si_1.append(vol_id_tot_old)
						moth_id_tot_Si_1.append(moth_id_tot_old)
						Strip_id_x_tot_Si_1.append(Strip_id_x_tot_old)
						Strip_id_y_tot_Si_1.append(Strip_id_y_tot_old)
						tray_id_tot_Si_1.append(tray_id_tot_old)
						plane_id_tot_Si_1.append(plane_id_tot_old)
						energy_dep_tot_Si_1.append(energy_dep_tot_old)					
					

						if where_other_vol != []:
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

			w_eth = where(energy_dep_tot >= E_th)
			where_eth = w_eth[0]


			event_id_tot = event_id_tot[where_eth]
			vol_id_tot = vol_id_tot[where_eth]
			moth_id_tot = moth_id_tot[where_eth]
			Strip_id_tot = Strip_id_tot[where_eth]
			Si_id_tot = Si_id_tot[where_eth]
			tray_id_tot = tray_id_tot[where_eth]
			plane_id_tot = plane_id_tot[where_eth]
			energy_dep_tot = energy_dep_tot[where_eth]
			pair_flag_tot = pair_flag_tot[where_eth]

			N_trig = len(np.uniqe(event_id_tot))
			event_array = event_id_tot[np.unique(event_id_tot)]


	







				# ------------------------------------
				
	ifile = ifile + 1



	













