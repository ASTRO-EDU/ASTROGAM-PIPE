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

			event_id_c = []
			track_id_c = []
			vol_id_c = []
			vol_name_c = []
			energy_dep_c = []
			x_en_c = []
			y_en_c = []
			z_en_c = []
			x_ex_c = []
			y_ex_c = []
			z_ex_c = []
			e_kin_en_c = []
			e_kin_ex_c = []
			en_dep_c = []
			
			theta_ent_c = []
			phi_ent_c = []
		
			theta_exit_c = []
			phi_exit_c = []

			child_id_c = []
			proc_id_c = []			
			moth_id_c = []

			if vol_id >= cal_vol_start and vol_id <= cal_vol_end:
                #print('b', energy_dep)
				if part_type == 'g' or energy_dep > 0.:
					event_id_calor = evt_id[i]
					moth_id_calor = mother_id[i]
					track_id_calor = trk_id[i]
					vol_id_calor = volume_id[i]
					vol_name_calor = volume_name[i]
					en_dep_calor = e_dep[i]
					x_en_calor = (x_ent[i])/10.
					y_en_calor = (y_ent[i])/10.
					z_en_calor = (z_ent[i])/10.
					x_ex_calor = (x_exit[i])/10.
					y_ex_calor = (y_exit[i])/10.
					z_ex_calor = (z_exit[i])/10.
					e_kin_en_calor = e_kin_ent[i]
					e_kin_ex_calor = e_kin_exit[i]
					en_dep_calor = e_dep[i]
				
					theta_ent_calor = (180./math.pi)*math.acos(-(cos_x_angle_ent))
					phi_ent_calor = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
		
					theta_exit_calor = (180./math.pi)*math.acos(-(cos_z_angle_exit))
					phi_exit_calor = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

					child_id_calor = parent_trk_id[i]
					proc_id_calor = process_id[i]

					event_id_c.append(event_id_calor)
					track_id_c.append(track_id_calor)
					vol_id_c.append(vol_id_calor)
					vol_name_c.append(vol_name_calor)
					energy_dep_c.append(en_dep_calor)
					x_en_c.append(x_en_calor)
					y_en_c.append(y_en_calor)
					z_en_c.append(z_en_calor)
					x_ex_c.append(x_ex_calor)
					y_ex_c.append(y_ex_calor)
					z_ex_c.append(z_ex_calor)
					e_kin_en_c.append(e_kin_en_calor)
					e_kin_ex_c.append(e_kin_ex_calor)
					en_dep_c.append(en_dep_calor)
				
					theta_ent_c.append(theta_ent_calor)
					phi_ent_c.append(phi_ent_calor)
		
					theta_exit_c.append(theta_exit_calor)
					phi_exit_c.append(phi_exit_calor)

					child_id_c.append(child_id_calor)
					proc_id_c.append(proc_id_calor)
					moth_id_c.append(moth_id_calor)

		# Reading the AC
		if ac_flag == 1:

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

			if vol_id >= ac_vol_start and vol_id <= ac_vol_end:
				if part_type == 'g' or energy_dep > 0.:
					event_id_anc = evt_id[i]
					vol_id_anc = vol_id[i]

					energy_dep_anc = e_dep[i]

					ent_x_anc = (x_ent[i])/10.
					ent_y_anc = (y_ent[i])/10.
					ent_z_anc = (z_ent[i])/10.
					exit_x_anc = (x_exit[i])/10.
					exit_y_anc = (y_exit[i])/10.
					exit_z_anc = (z_exit[i])/10.

					theta_ent_anc = (180./math.pi)*math.acos(-(cos_x_angle_ent))
					phi_ent_anc = (180./math.pi)*math.atan((cos_y_angle_ent)/(cos_x_angle_ent))
	
					theta_exit_anc = (180./math.pi)*math.acos(-(cos_z_angle_exit))
					phi_exit_anc = (180./math.pi)*math.atan((cos_y_angle_exit)/(cos_x_angle_exit))

					child_id_anc = parent_trk_id[i]
					proc_id_anc = process_id[i]

					if isStrip == 1:
						moth_id_anc = mother_id[i]
					else:
						moth_id_anc = 0

					event_id_ac.append(event_id_anc)
					moth_id_ac.append(moth_id_anc)
					vol_id_ac.append(vol_id_anc)
					x_en_ac.append(x_en_anc)
					y_en_ac.append(y_en_anc)
					z_en_ac.append(z_en_anc)
					x_ex_ac.append(x_ex_anc)
					y_ex_ac.append(y_ex_anc)
					z_ex_ac.append(z_ex_anc)
				
					theta_ent_ac.append(theta_ent_anc)
					phi_ent_ac.append(phi_ent_anc)
		
					theta_exit_ac.append(theta_exit_anc)
					phi_exit_ac.append(phi_exit_anc)

					child_id_ac.append(child_id_anc)
					proc_id_ac.append(proc_id_anc)

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

					if len(where_pair) != 0:
						pair_flag_tot_old = 1
					else:
						pair_flag_tot_old = 0				#il terzultimo valore di IDL e' 1 quando dovrebbe essere 0 perche nel where_pair e' -1   
					pair_flag_tot.append(pair_flag_tot_old)	
 
					w_compton = np.where((all_child == 1) & (all_proc == 3))
					where_compton = w_compton[0]

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
				#j_max = max(where_event_eq)
				#j = j_max + 1
			
				N_event_eq = len(where_event_eq)                                          
				if where_event_eq[N_event_eq-1] < len(event_id_tr)-1:
					j = where_event_eq[N_event_eq-1]+1					
				else:
					break

			event_id_tot_temp = [0]*(2*len(event_id_tot))
			vol_id_tot_temp = [0]*(2*len(event_id_tot))
			moth_id_tot_temp = [0]*((2*len(event_id_tot)))
			Strip_id_tot_temp = [0]*(2*len(event_id_tot))
			Si_id_tot_temp = [0]*((2*len(event_id_tot)))
			tray_id_tot_temp = [0]*(2*len(event_id_tot))
			plane_id_tot_temp = [0]*(2*len(event_id_tot))
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

			e_dep_temp_list = [] 		
			event_id_tot_list = []
			vol_id_tot_list = []
			moth_id_tot_list = []
			Strip_id_tot_list = []
			Si_id_tot_list = []
			tray_id_tot_list = []
			plane_id_tot_list = []
			energy_dep_tot_list = []
			pair_flag_tot_list = []


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

						e_dep_temp_list.append(e_dep_temp_old)
						event_id_tot_list.append(event_id_tot_old)
						vol_id_tot_list.append(vol_id_tot_old)
						moth_id_tot_list.append(moth_id_tot_old)
						Strip_id_tot_list.append(Strip_id_tot_old)
						Si_id_tot_list.append(Si_id_tot_old)
						tray_id_tot_list.append(tray_id_tot_old)
						plane_id_tot_list.append(plane_id_tot_old)
						energy_dep_tot_list.append(energy_dep_tot_old)					
						pair_flag_tot_list.append(pair_flag_tot_old)

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


					vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]) & (Si_id_temp == 1))	
					where_vol_eq = vol_eq[0]
					
					other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]) | (Si_id_temp != 1))
					where_other_vol = other_vol[0]
					
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

						e_dep_temp_list.append(e_dep_temp_old)
						event_id_tot_list.append(event_id_tot_old)
						vol_id_tot_list.append(vol_id_tot_old)
						moth_id_tot_list.append(moth_id_tot_old)
						Strip_id_tot_list.append(Strip_id_tot_old)
						Si_id_tot_list.append(Si_id_tot_old)
						tray_id_tot_list.append(tray_id_tot_old)
						plane_id_tot_list.append(plane_id_tot_old)
						energy_dep_tot_list.append(energy_dep_tot_old)					
						pair_flag_tot_list.append(pair_flag_tot_old)

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
			e_dep_temp = np.array(e_dep_temp_list) 		
			event_id_tot = np.array(event_id_tot_list)
			vol_id_tot = np.array(vol_id_tot_list)
			moth_id_tot = np.array(moth_id_tot_list)
			Strip_id_tot = np.array(Strip_id_tot_list)
			Si_id_tot = np.array(Si_id_tot_list)
			tray_id_tot = np.array(tray_id_tot_list)
			plane_id_tot = np.array(plane_id_tot_list)
			energy_dep_tot = np.array(energy_dep_tot_list)
			pair_flag_tot = np.array(pair_flag_tot_list)

			w_eth = np.where(energy_dep_tot >= E_th)
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


			####Funzione indici unici event_id_tot

			event_id_index = []

			for i in range(len(event_id_tot)):
				event_id_index_old = i
				event_id_index.append(i)

			event_uniq_index = []	

			i = 1
			while i < len(event_id_tot):

				if event_id_tot[i] != event_id_tot[i-1]:
					event_uniq_index.append(event_id_index[i-1])

				if i == len(event_id_tot)-1:
					event_uniq_index.append(event_id_index[i])

				i = i + 1


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
				w_event_eq = np.where(event_id_tot == event_id_tot[j])
				where_event_eq = w_event_eq[0]

				event_id_temp = event_id_tot[where_event_eq]
				vol_id_temp = vol_id_tot[where_event_eq]
				moth_id_temp = moth_id_tot[where_event_eq]
				Strip_id_temp = Strip_id_tot[where_event_eq]
				Si_id_temp = Si_id_tot[where_event_eq]
				tray_id_temp = tray_id_tot[where_event_eq]
				plane_id_temp = plane_id_tot[where_event_eq]
				energy_dep_temp = energy_dep_tot[where_event_eq]
				pair_flag_temp = pair_flag_tot[where_event_eq]			
				
				# Funzione di ordinamento

				vol_sort_arr_list = []
				vol_id_temp_sort = vol_id_temp
				for z in range(len(vol_id_temp_sort)):

					vol_sort_arr_old = z
					vol_sort_arr_list.append(vol_sort_arr_old)
						
				swap = True
				while swap == True:
					i = 1
					swap = False
												
					while i < len(vol_id_temp_sort):
						if vol_id_temp_sort[i-1] > vol_id_temp_sort[i]:
							temp = vol_id_temp_sort[i-1]
							vol_id_temp_sort[i-1] = vol_id_temp_sort[i]
							vol_id_temp_sort[i] = temp
							temp = vol_sort_arr_list[i-1]
							vol_sort_arr_list[i-1] = vol_sort_arr_list[i]
							vol_sort_arr_list[i] = temp
							swap = True
						i = i+1
						
				vol_sort_arr = np.array(vol_sort_arr_list)				
				
				# fine funzione di ordinamento	

				vol_id_temp = np.ravel(vol_id_temp)
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


			Glob_event_id_test_list = []
			Glob_vol_id_test_list = []
			Glob_moth_id_test_list = []
			Glob_Strip_id_test_list = []
			Glob_Si_id_test_list = []
			Glob_tray_id_test_list = []
			Glob_plane_id_test_list = []
			Glob_pos_test_list = []
			Glob_zpos_test_list = []
			Glob_energy_dep_test_list = []
			Glob_pair_flag_test_list = []


			for j in range(N_trig):

				where_test_x_old = np.where(Glob_energy_dep_x_top[:,j] > 0.)
				where_test_x = where_test_x_old[0]
				

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

				where_test_y_old = np.where(Glob_energy_dep_y_top[:,j] > 0.)
				where_test_y = where_test_y_old[0]

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
										

				
				
				# Funzione di ordinamento

				tray_sort_arr_list = []
				Glob_tray_id_test_temp_sort = Glob_tray_id_test_temp
				for z in range(len(Glob_tray_id_test_temp_sort)):

					tray_sort_arr_old = z
					tray_sort_arr_list.append(tray_sort_arr_old)
						
				swap = True
				while swap == True:
					i = 1
					swap = False
												
					while i < len(Glob_tray_id_test_temp_sort):
						if Glob_tray_id_test_temp_sort[i-1] > Glob_tray_id_test_temp_sort[i]:
							temp = Glob_tray_id_test_temp_sort[i-1]
							Glob_tray_id_test_temp_sort[i-1] = Glob_tray_id_test_temp_sort[i]
							Glob_tray_id_test_temp_sort[i] = temp
							temp = tray_sort_arr_list[i-1]
							tray_sort_arr_list[i-1] = tray_sort_arr_list[i]
							tray_sort_arr_list[i] = temp
							swap = True
						i = i+1
						
				tray_sort_arr_temp = np.array(tray_sort_arr_list)				

				# fine funzione di ordinamento	

				tray_sort_arr = tray_sort_arr_temp[::-1]

				Glob_vol_id_test_temp = Glob_vol_id_test_temp[tray_sort_arr]
			        Glob_moth_id_test_temp = Glob_moth_id_test_temp[tray_sort_arr]
				Glob_Strip_id_test_temp = Glob_Strip_id_test_temp[tray_sort_arr]
				Glob_Si_id_test_temp = Glob_Si_id_test_temp[tray_sort_arr]
				Glob_tray_id_test_temp_descending = sorted(Glob_tray_id_test_temp[tray_sort_arr], reverse=True)
				Glob_tray_id_test_temp = np.array(Glob_tray_id_test_temp_descending)
				Glob_plane_id_test_temp = Glob_plane_id_test_temp[tray_sort_arr]
				Glob_pos_test_temp = Glob_pos_test_temp[tray_sort_arr]
				Glob_zpos_test_temp = Glob_zpos_test_temp[tray_sort_arr]
				Glob_energy_dep_test_temp = Glob_energy_dep_test_temp[tray_sort_arr]
				Glob_pair_flag_test_temp = Glob_pair_flag_test_temp[tray_sort_arr]


				vol_id_intray_list = []
				moth_id_intray_list = []
				Strip_id_intray_list = []
				Si_id_intray_list = []
				tray_id_intray_list = []
				plane_id_intray_list = []
				pos_intray_list = []
				zpos_intray_list = []
				energy_dep_intray_list = []
				pair_flag_intray_list = []


		
				intray = 0
				while 1:
					where_tray_eq_old = np.where(Glob_tray_id_test_temp == Glob_tray_id_test_temp[intray])
					where_tray_eq = where_tray_eq_old[0]

					where_other_tray_old = np.where(Glob_tray_id_test_temp != Glob_tray_id_test_temp[intray])
					where_other_tray = where_other_tray_old[0]


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

					where_Y_old = np.where(Si_id_extract == 1)
					where_Y = where_Y_old[0] 

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

						vol_id_intray_list.append(vol_id_intray_old)
						moth_id_intray_list.append(moth_id_intray_old)
						Strip_id_intray_list.append(Strip_id_intray_old)
						Si_id_intray_list.append(Si_id_intray_old)
						tray_id_intray_list.append(tray_id_intray_old)
						plane_id_intray_list.append(plane_id_intray_old)
						pos_intray_list.append(pos_intray_old)
						zpos_intray_list.append(zpos_intray_old)
						energy_dep_intray_list.append(energy_dep_intray_old)
						pair_flag_intray_list.append(pair_flag_intray_old)


					where_X_old = np.where(Si_id_extract == 0)
					where_X = where_X_old[0]

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

						vol_id_intray_list.append(vol_id_intray_old)
						moth_id_intray_list.append(moth_id_intray_old)
						Strip_id_intray_list.append(Strip_id_intray_old)
						Si_id_intray_list.append(Si_id_intray_old)
						tray_id_intray_list.append(tray_id_intray_old)
						plane_id_intray_list.append(plane_id_intray_old)
						pos_intray_list.append(pos_intray_old)
						zpos_intray_list.append(zpos_intray_old)
						energy_dep_intray_list.append(energy_dep_intray_old)
						pair_flag_intray_list.append(pair_flag_intray_old)

 
					N_tray_eq = len(where_tray_eq)
					if where_tray_eq[N_tray_eq-1] < len(Glob_tray_id_test_temp)-1:
						intray = where_tray_eq[N_tray_eq-1]+1						
					else:
						break 
 
				vol_id_intray_concatenate = np.ma.concatenate(vol_id_intray_list)
				moth_id_intray_concatenate = np.ma.concatenate(moth_id_intray_list)
				Strip_id_intray_concatenate = np.ma.concatenate(Strip_id_intray_list)
				Si_id_intray_concatenate = np.ma.concatenate(Si_id_intray_list)
				tray_id_intray_concatenate = np.ma.concatenate(tray_id_intray_list)
				plane_id_intray_concatenate = np.ma.concatenate(plane_id_intray_list)
				pos_intray_concatenate = np.ma.concatenate(pos_intray_list)
				zpos_intray_concatenate = np.ma.concatenate(zpos_intray_list)
				energy_dep_intray_concatenate = np.ma.concatenate(energy_dep_intray_list)
				pair_flag_intray_concatenate = np.ma.concatenate(pair_flag_intray_list)
 				
				vol_id_temp = np.array(vol_id_intray_concatenate)
				moth_id_temp = np.array(moth_id_intray_concatenate)
				Strip_id_temp = np.array(Strip_id_intray_concatenate)
				Si_id_temp = np.array(Si_id_intray_concatenate)
				tray_id_temp = np.array(tray_id_intray_concatenate)
				plane_id_temp = np.array(plane_id_intray_concatenate)
				pos_temp = np.array(pos_intray_concatenate)
				zpos_temp = np.array(zpos_intray_concatenate)
				energy_dep_temp = np.array(energy_dep_intray_concatenate)
				pair_flag_temp = np.array(pair_flag_intray_concatenate)				  
 
 
				event_id_temp = []
				for k in range(len(vol_id_temp)):
					event_id_temp_old = event_array[j]
					event_id_temp.append(event_id_temp_old)



				Glob_event_id_test_old = event_id_temp
				Glob_vol_id_test_old = vol_id_temp
				Glob_moth_id_test_old = moth_id_temp
				Glob_Strip_id_test_old = Strip_id_temp
				Glob_Si_id_test_old = Si_id_temp
				Glob_tray_id_test_old = tray_id_temp
				Glob_plane_id_test_old = plane_id_temp
				Glob_pos_test_old = pos_temp
				Glob_zpos_test_old = zpos_temp
				Glob_energy_dep_test_old = energy_dep_temp
				Glob_pair_flag_test_old = pair_flag_temp


				Glob_event_id_test_list.append(Glob_event_id_test_old)
				Glob_vol_id_test_list.append(Glob_vol_id_test_old)
				Glob_moth_id_test_list.append(Glob_moth_id_test_old)
				Glob_Strip_id_test_list.append(Glob_Strip_id_test_old)
				Glob_Si_id_test_list.append(Glob_Si_id_test_old)
				Glob_tray_id_test_list.append(Glob_tray_id_test_old)
				Glob_plane_id_test_list.append(Glob_plane_id_test_old)
				Glob_pos_test_list.append(Glob_pos_test_old)
				Glob_zpos_test_list.append(Glob_zpos_test_old)
				Glob_energy_dep_test_list.append(Glob_energy_dep_test_old)
				Glob_pair_flag_test_list.append(Glob_pair_flag_test_old)
	
			Glob_event_id_test_concatenate = np.ma.concatenate(Glob_event_id_test_list)	
			Glob_vol_id_test_concatenate = np.ma.concatenate(Glob_vol_id_test_list)
			Glob_moth_id_test_concatenate = np.ma.concatenate(Glob_moth_id_test_list)
			Glob_Strip_id_test_concatenate = np.ma.concatenate(Glob_Strip_id_test_list)
			Glob_Si_id_test_concatenate = np.ma.concatenate(Glob_Si_id_test_list)
			Glob_tray_id_test_concatenate = np.ma.concatenate(Glob_tray_id_test_list)
			Glob_plane_id_test_concatenate = np.ma.concatenate(Glob_plane_id_test_list)
			Glob_pos_test_concatenate = np.ma.concatenate(Glob_pos_test_list)
			Glob_zpos_test_concatenate = np.ma.concatenate(Glob_zpos_test_list)
			Glob_energy_dep_test_concatenate = np.ma.concatenate(Glob_energy_dep_test_list)
			Glob_pair_flag_test_concatenate = np.ma.concatenate(Glob_pair_flag_test_list)			

			Glob_event_id_test = np.array(Glob_event_id_test_concatenate)
			Glob_vol_id_test = np.array(Glob_vol_id_test_concatenate)
			Glob_moth_id_test = np.array(Glob_moth_id_test_concatenate)
			Glob_Strip_id_test = np.array(Glob_Strip_id_test_concatenate)
			Glob_Si_id_test = np.array(Glob_Si_id_test_concatenate)
			Glob_tray_id_test = np.array(Glob_tray_id_test_concatenate)
			Glob_plane_id_test = np.array(Glob_plane_id_test_concatenate)
			Glob_pos_test = np.array(Glob_pos_test_concatenate)
			Glob_zpos_test = np.array(Glob_zpos_test_concatenate)
			Glob_energy_dep_test = np.array(Glob_energy_dep_test_concatenate)
			Glob_pair_flag_test = np.array(Glob_pair_flag_test_concatenate, dtype = np.int64)

			#print(Glob_event_id_test_list)
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





			#print(Glob_event_id_test)

			j=0
			while 1:
				where_event_eq_arr = np.where(Glob_event_id_test == Glob_event_id_test[j])
				where_event_eq = where_event_eq_arr[0]
				print(j)

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

				where_x_arr = np.where(Glob_Si_id_test_temp == 0)
				where_x = where_x_arr[0]				
				
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

				where_y_arr = np.where(Glob_Si_id_test_temp == 1)
				where_y = where_y_arr[0]				
				
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


			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('       L0.5 - cluster baricenter ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


				# ------------------------------------
				
	ifile = ifile + 1



	













