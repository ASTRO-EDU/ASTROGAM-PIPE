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

class DHSim:
### parameter function ###########

    # PARAMETERS DEFINITION: it keeps the input parameters from eASTROGAM_ANALYSIS_INPUT AND RETURNS THE STRINGS FOR THE OUTPUTS #

    def parameter(self, sim_type, py_list, part_type, ene_range, pol_type, source_g, isStrip, repli, cal_flag, ac_flag, astrogam_version, bogemms_tag, ene_min, ene_max , passive_flag, energy_thresh):

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

        if astrogam_version=='V1.0':
			astrogam_tag = '01'
			sim_tag = 'eAST'+str(bogemms_tag)+str(astrogam_tag)+'0102'

        if astrogam_version=='V1.1':
			astrogam_tag = '11'
			sim_tag = 'eAST'+str(bogemms_tag)+str(astrogam_tag)+'2021'

        if astrogam_version=='V1.2':
			astrogam_tag = '12'
			sim_tag = 'eAST'+str(bogemms_tag)+str(astrogam_tag)+'2022'

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

        if astrogam_version == 'V1.0' or astrogam_version == 'V1.1' or astrogam_version == 'V1.2' or astrogam_version == 'V2.0' or astrogam_version=='V10.0':
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



		# setting specific agile version variables
        if astrogam_version=='V1.0' or astrogam_version=='V1.1' or astrogam_version=='V1.2' or astrogam_version=='V2.0' or astrogam_version=='V10.0':
            # --------> volume ID
            tracker_top_vol_start = 1090000
            tracker_bottom_vol_start = 1000000
            tracker_top_bot_diff = 90000
            # for the S1 format
            x_layer_id = 90000
            y_layer_id = 80000


            if astrogam_version=='V1.0':
                cal_vol_start = 50000
                cal_vol_end = 58463
                # --------> design
                N_tray = 56
                N_plane = N_tray*1
                N_strip = 3840
                tray_side = 92.16 #cm
                strip_side = tray_side/N_strip

            if astrogam_version=='V1.1':
                cal_vol_start = 50000
                cal_vol_end = 83855
                N_tray = 56
                # --------> design
                N_tray = 56
                N_plane = N_tray*1
                N_strip = 3840
                tray_side = 92.16 #cm
                strip_side = tray_side/N_strip

            if astrogam_version=='V1.2':
                cal_vol_start = 50000
                cal_vol_end = 53135
                # --------> design
                N_tray = 65
                N_plane = N_tray*1
                N_strip = 1152
                tray_side = 55.296 #cm
                strip_side = tray_side/N_strip

            if astrogam_version=='V2.0':
                cal_vol_start = 50000
                cal_vol_end = 50783
                # --------> design
                N_tray = 25
                N_plane = N_tray*1
                N_strip = 1152
                tray_side = 27.648 #cm
                strip_side = tray_side/N_strip

            if astrogam_version=='V10.0':
                cal_vol_start = 50000
                cal_vol_end = 51443
                # --------> design
                N_tray = 60
                N_plane = N_tray*1
                N_strip = 760
                tray_side = 38.0 #cm
                strip_side = tray_side/N_strip

            ac_vol_start = 301
            ac_vol_end = 350

            panel_S = [301, 302, 303]
            panel_D = [311, 312, 313]
            panel_F = [321, 322, 323]
            panel_B = [331, 332, 333]
            panel_top = 340


            # --------> processing
            # accoppiamento capacitivo
            # acap = [0.035, 0.045, 0.095, 0.115, 0.38, 1., 0.38, 0.115, 0.095, 0.045, 0.035]

            # tracker energy threshold (0.25 MIP)
            E_th = float(energy_thresh)  # keV

            E_th_cal = 30. # keV
		
		
        self.astrogam_tag = astrogam_tag
        self.sim_tag = sim_tag
        self.ene_type = ene_type
        self.py_dir = py_dir
        self.py_name = py_name
        self.sim_name = sim_name
        self.pol_string = pol_string
        self.sdir = sdir
        self.sname = sname
        self.dir_cal = dir_cal
        self.dir_passive = dir_passive
        self.stripDir = stripDir
        self.stripname = stripname
        self.tracker_top_vol_start = tracker_top_vol_start
        self.tracker_bottom_vol_start = tracker_bottom_vol_start
        self.tracker_top_bot_diff = tracker_top_bot_diff
        self.x_layer_id = x_layer_id
        self.y_layer_id = y_layer_id
        self.cal_vol_start = cal_vol_start
        self.cal_vol_end = cal_vol_end
        self.ac_vol_start = ac_vol_start
        self.ac_vol_end = ac_vol_end
        self.panel_S = panel_S
        self.panel_D = panel_D
        self.panel_F = panel_F
        self.panel_B = panel_B
        self.panel_top = panel_top
        self.N_tray = N_tray
        self.N_plane = N_plane
        self.N_strip = N_strip
        self.tray_side = tray_side
        self.strip_side = strip_side
        self.E_th = E_th
        self.E_th_cal = E_th_cal
        self.ene_dis = ene_dis
 		
####### FILEPATH #########
    # IT RETURNS THE SIMULATION PATH #
    def get_filepath(self, astrogam_version, theta_type, N_in, part_type):
	
		filepath = './input_eASTROGAM'+astrogam_version+self.sdir+'/theta'+str(theta_type)+'/'+self.stripDir+self.py_dir+self.dir_cal+self.dir_passive+'/'+self.ene_type+'MeV.'+self.sim_name+'.'+str(theta_type)+'theta.'+self.pol_string+str(N_in)+part_type

		self.filepath = filepath
	
		print('eASTROGAM simulation path:' + filepath)


########## OUTDIR ######
    # IT RETURNS THE PATH OF THE OUTPUT DIRECTORY #
    def get_outdir(self, astrogam_version, theta_type, N_in, part_type, energy_thresh):


		outdir = ('./output_eASTROGAM'+astrogam_version+self.sdir+'/theta'+str(theta_type)+'/'+self.stripDir+self.py_dir+'/'+str(self.sim_name)+'/'+str(self.ene_type)+'MeV/'+str(N_in)+part_type+self.dir_cal+self.dir_passive+'/'+str(energy_thresh)+'keV')

		self.outdir = outdir

		if not os.path.exists(outdir):
			out_dir = os.makedirs(outdir,0777)
		
##### READING FITS ###########
    # IT READS THE FITS FILE GENERATED BY BOGEMMS AND CREATES THE ARRAYS FOR THE TRACKER, THE CALORIMETER AND THE ANTICOINCIDENCE CONSIDERING ONLY THE EVENTS WITH DEPOSITED ENERGY > 0#
    def reading_fits(self, ifile, cal_flag, ac_flag, part_type, isStrip):

		t = fits.open(self.filepath+'/xyz.'+str(ifile)+'.fits.gz')
   
		tbdata = t[1].data

		evt_id = tbdata.field('EVT_ID')
		trk_id = tbdata.field('TRK_ID')
		parent_trk_id = tbdata.field('PARENT_TRK_ID')
		volume_id = tbdata.field('VOLUME_ID')
		mother_id = tbdata.field('MOTHER_ID')
		e_dep = tbdata.field('E_DEP')
		x_ent = tbdata.field('X_ENT')
		y_ent = tbdata.field('Y_ENT')
		z_ent = tbdata.field('Z_ENT')
		x_exit = tbdata.field('X_EXIT')
		y_exit = tbdata.field('Y_EXIT')
		z_exit = tbdata.field('Z_EXIT')
		mdx_ent = tbdata.field('MDX_ENT')
		mdy_ent = tbdata.field('MDY_ENT')
		mdz_ent = tbdata.field('MDZ_ENT')
		mdx_exit = tbdata.field('MDX_EXIT')
		mdy_exit = tbdata.field('MDY_EXIT')
		mdz_exit = tbdata.field('MDZ_EXIT')
		gtime_ent = tbdata.field('GTIME_ENT')
		particle_id = tbdata.field('PARTICLE_ID')
		process_id = tbdata.field('PROCESS_ID')
	
	
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
		part_id_tr = []
		trk_id_tr = []
		child_id_tr = []
		proc_id_tr = []
		gtime_ent_tr = []
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
		part_id_cal = []
		trk_id_cal = []
		child_id_cal = []
		proc_id_cal = []
		gtime_ent_cal = []
	
		event_id_ac = []
		vol_id_ac = []
		energy_dep_ac = []
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
		part_id_ac = []
		trk_id_ac = []
		child_id_ac = []
		proc_id_ac = []
		gtime_ent_ac = []
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
		
			#  Reading the tracker (events with E > 0)
	       		        
			if vol_id >= self.tracker_bottom_vol_start or moth_id >= self.tracker_bottom_vol_start:

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
					part_id_tr.append(particle_id[i])
					trk_id_tr.append(trk_id[i])
					gtime_ent_tr.append(gtime_ent[i])					


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
					part_id_tr.append(particle_id[i])
					trk_id_tr.append(trk_id[i])
					gtime_ent_tr.append(gtime_ent[i])					
				

			# Reading the Calorimeter
			if cal_flag == 1:

				if vol_id >= self.cal_vol_start and vol_id <= self.cal_vol_end:
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
						part_id_cal.append(particle_id[i])
						trk_id_cal.append(trk_id[i])
						gtime_ent_cal.append(gtime_ent[i])					
					

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
						part_id_cal.append(particle_id[i])
						trk_id_cal.append(trk_id[i])
						gtime_ent_cal.append(gtime_ent[i])					


			# Reading the AC
			if ac_flag == 1:

				if vol_id >= self.ac_vol_start and vol_id <= self.ac_vol_end:
	
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
						vol_id_ac.append(volume_id[i])
						x_en_ac.append(x_ent[i]/10.)
						y_en_ac.append(y_ent[i]/10.)
						z_en_ac.append(z_ent[i]/10.)
						x_ex_ac.append(x_exit[i]/10.)
						y_ex_ac.append(y_exit[i]/10.)
						z_ex_ac.append(z_exit[i]/10.)
						energy_dep_ac.append(e_dep[i])
						theta_ent_ac.append(theta_ent_anc)
						phi_ent_ac.append(phi_ent_anc)
		
						theta_exit_ac.append(theta_exit_anc)
						phi_exit_ac_.append(phi_exit_anc)

						child_id_ac.append(parent_trk_id[i])
						proc_id_ac.append(process_id[i])
						part_id_ac.append(particle_id[i])
						trk_id_ac.append(trk_id[i])
						gtime_ent_ac.append(gtime_ent[i])					


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
						vol_id_ac.append(volume_id[i])
						x_en_ac.append(x_ent[i]/10.)
						y_en_ac.append(y_ent[i]/10.)
						z_en_ac.append(z_ent[i]/10.)
						x_ex_ac.append(x_exit[i]/10.)
						y_ex_ac.append(y_exit[i]/10.)
						z_ex_ac.append(z_exit[i]/10.)
						energy_dep_ac.append(e_dep[i])
						theta_ent_ac.append(theta_ent_anc)
						phi_ent_ac.append(phi_ent_anc)
		
						theta_exit_ac.append(theta_exit_anc)
						phi_exit_ac.append(phi_exit_anc)

						child_id_ac.append(parent_trk_id[i])
						proc_id_ac.append(process_id[i])
						part_id_ac.append(particle_id[i])
						trk_id_ac.append(trk_id[i])
						gtime_ent_ac.append(gtime_ent[i])					

	

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
		theta_ent_tr = np.array(theta_ent_tr)
		phi_ent_tr = np.array(phi_ent_tr)
		theta_exit_tr = np.array(theta_exit_tr)
		phi_exit_tr = np.array(phi_exit_tr)
		child_id_tr = np.array(child_id_tr)
		proc_id_tr = np.array(proc_id_tr)
		x_pos = np.array(x_pos)
		y_pos = np.array(y_pos)
		z_pos = np.array(z_pos)
		part_id_tr = np.array(part_id_tr)
		trk_id_tr = np.array(trk_id_tr)
		gtime_ent_tr = np.array(gtime_ent_tr)					

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
		part_id_cal = np.array(part_id_cal)
		trk_id_cal = np.array(trk_id_cal)
		gtime_ent_cal = np.array(gtime_ent_cal)					


		#### AC arrays

		event_id_ac = np.array(event_id_ac)
		vol_id_ac = np.array(vol_id_ac)
		energy_dep_ac = np.array(energy_dep_ac)
		x_en_ac = np.array(x_en_ac)
		y_en_ac = np.array(y_en_ac)
		z_en_ac = np.array(z_en_ac)
		x_ex_ac = np.array(x_ex_ac)
		y_ex_ac = np.array(y_ex_ac)
		z_ex_ac = np.array(z_ex_ac)
		theta_ent_ac = np.array(theta_ent_ac)
		phi_ent_ac = np.array(phi_ent_ac)
		theta_exit_ac = np.array(theta_exit_ac)
		phi_exit_ac = np.array(phi_exit_ac)
		child_id_ac = np.array(child_id_ac)
		proc_id_ac = np.array(proc_id_ac)			
		moth_id_ac = np.array(moth_id_ac)
		part_id_ac = np.array(part_id_ac)
		trk_id_ac = np.array(trk_id_ac)
		gtime_ent_ac = np.array(gtime_ent_ac)					

		self.vol_id_tr = vol_id_tr 
		self.moth_id_tr = moth_id_tr 
		self.event_id_tr = event_id_tr 
		self.en_dep_tr = en_dep_tr 
		self.x_en_tr = x_en_tr 
		self.y_en_tr = y_en_tr
		self.z_en_tr = z_en_tr 
		self.x_ex_tr = x_ex_tr 
		self.y_ex_tr = y_ex_tr 
		self.z_ex_tr = z_ex_tr
		self.theta_ent_tr = theta_ent_tr
		self.phi_ent_tr = phi_ent_tr
		self.theta_exit_tr = theta_exit_tr
		self.phi_exit_tr = phi_exit_tr
		self.child_id_tr = child_id_tr 
		self.proc_id_tr = proc_id_tr 
		self.theta_ent_tr = theta_ent_tr 
		self.phi_ent_tr = phi_ent_tr 
		self.theta_exit_tr = theta_exit_tr 
		self.phi_exit_tr = phi_exit_tr 
		self.x_pos = x_pos 
		self.y_pos = y_pos
		self.z_pos = z_pos 
		self.part_id_tr = part_id_tr
		self.trk_id_tr = trk_id_tr
		self.gtime_ent_tr = gtime_ent_tr					

		self.vol_id_cal = vol_id_cal 
		self.moth_id_cal = moth_id_cal
		self.event_id_cal = event_id_cal
		self.energy_dep_cal = energy_dep_cal
		self.x_en_cal = x_en_cal
		self.y_en_cal = y_en_cal
		self.z_en_cal = z_en_cal
		self.x_ex_cal = x_ex_cal 
		self.y_ex_cal = y_ex_cal 
		self.z_ex_cal = z_ex_cal 
		self.child_id_cal = child_id_cal 
		self.proc_id_cal = proc_id_cal 
		self.theta_ent_cal = theta_ent_cal 
		self.phi_ent_cal = phi_ent_cal 
		self.theta_exit_cal = theta_exit_cal 
		self.phi_exit_cal = phi_exit_cal 
		self.part_id_cal = part_id_cal
		self.trk_id_cal = trk_id_cal
		self.gtime_ent_cal = gtime_ent_cal					

		self.vol_id_ac = vol_id_ac 
		self.event_id_ac = event_id_ac
		self.energy_dep_ac = energy_dep_ac 
		self.x_en_ac = x_en_ac 
		self.y_en_ac = y_en_ac 
		self.z_en_ac = z_en_ac 
		self.x_ex_ac = x_ex_ac 
		self.y_ex_ac = y_ex_ac 
		self.z_ex_ac = z_ex_ac 
		self.child_id_ac = child_id_ac 
		self.proc_id_ac = proc_id_ac 
		self.theta_ent_ac = theta_ent_ac 
		self.phi_ent_ac = phi_ent_ac 
		self.theta_exit_ac = theta_exit_ac 
		self.phi_exit_ac = phi_exit_ac
		self.part_id_ac = part_id_ac
		self.trk_id_ac = trk_id_ac
		self.gtime_ent_ac = gtime_ent_ac					
		self.moth_id_ac = moth_id_ac
		

########## CONVERSION TRAY AND PLANE
    # IT CREATES THE STRIP ID, PLANE ID AND TRAY ID ARRAYS #
    def conversion(self, isStrip, repli):

		# From Tracker volume ID to strip

		Strip_id_x = []
		Strip_id_y = []
		tray_id = []
		
		# Conversion from tray ID (starting from bottom) to plane ID (starting from the top)
		plane_id = []

		j=0
		
		while j < len(self.vol_id_tr):
		
			if isStrip == 1:
				if repli == 1:
					Strip_y = self.vol_id_tr[j]
					tray = self.moth_id_tr[j]/self.tracker_bottom_vol_start
					invert_tray_id = (self.N_tray - tray)+1
					vol_id_temp = self.moth_id_tr[j] - (self.tracker_bottom_vol_start*tray + self.tracker_top_bot_diff) # removing 1000000xn_tray + 90000
					Strip_x = vol_id_temp
					plane = invert_tray_id	

					plane_id.append(plane)
					Strip_id_y.append(Strip_y)
					Strip_id_x.append(Strip_x)
					tray_id.append(tray)

			else:	
				
				Strip_y = 0
				tray = self.vol_id_tr[j]/self.tracker_bottom_vol_start
				invert_tray_id = (self.N_tray - tray)+1
				Strip_x= 0
				plane = invert_tray_id				
				
				plane_id.append(plane)
				Strip_id_y.append(Strip_y)
				Strip_id_x.append(Strip_x)
				tray_id.append(tray)

			j = j+1

		tray_id = np.array(tray_id)
		plane_id = np.array(plane_id)
		Strip_id_x = np.array(Strip_id_x)
		Strip_id_y = np.array(Strip_id_y)

		self.Strip_id_x = Strip_id_x 
		self.Strip_id_y = Strip_id_y 
		self.tray_id = tray_id 
		self.plane_id = plane_id

	
####### WRITING #########
    ###### G4_raw.fits  ########
    # IT WRITES THE FITS FILE WITH RAW DATA OF THE TRACKER #
    
    def writing_G4raw(self, N_in, part_type, theta_type, phi_type, ifile, astrogam_version):

		
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tr)	
		col2 = fits.Column(name='VOL_ID', format='I', array=self.vol_id_tr)
		col3 = fits.Column(name='MOTH_ID', format='J', array=self.moth_id_tr)
		col4 = fits.Column(name='TRAY_ID', format='I', array=self.tray_id)
		col5 = fits.Column(name='PLANE_ID', format='I', array=self.plane_id)
		col6 = fits.Column(name='STRIP_ID_X', format='I', array=self.Strip_id_x)
		col7 = fits.Column(name='STRIP_ID_Y', format='I', array=self.Strip_id_y)
		col8 = fits.Column(name='E_DEP', format='F20.5', array=self.en_dep_tr)
		col9 = fits.Column(name='X_ENT', format='F20.5', array=self.x_en_tr)
		col10 = fits.Column(name='Y_ENT', format='F20.5', array=self.y_en_tr)
		col11 = fits.Column(name='Z_ENT', format='F20.5', array=self.z_en_tr)
		col12 = fits.Column(name='X_EXIT', format='F20.5', array=self.x_ex_tr)
		col13 = fits.Column(name='Y_EXIT', format='F20.5', array=self.y_ex_tr)
		col14 = fits.Column(name='Z_EXIT', format='F20.5', array=self.z_ex_tr)
		col15 = fits.Column(name='PART_ID', format='I', array=self.part_id_tr)
		col16 = fits.Column(name='TRK_ID', format='I', array=self.trk_id_tr)
		col17 = fits.Column(name='CHILD_ID', format='I', array=self.child_id_tr)
		col18 = fits.Column(name='PROC_ID', format='I', array=self.proc_id_tr)

		
		cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy ='+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Position unit = cm', ext=1)

		fits.setval(self.outdir+'/G4.RAW.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = keV', ext=1)
		
    #### L0.fits  ############
    # IT WRITES THE FITS FILE WITH STRIP DATA OF THE TRACKER #

    def writing_L0(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):


		col1 = fits.Column(name='EVT_ID', format='I', array=self.Glob_event_id_test)	
		col2 = fits.Column(name='VOL_ID', format='J', array=self.Glob_vol_id_test)
		col3 = fits.Column(name='MOTH_ID', format='J', array=self.Glob_moth_id_test)
		col4 = fits.Column(name='TRAY_ID', format='I', array=self.Glob_tray_id_test)
		col5 = fits.Column(name='PLANE_ID', format='I', array=self.Glob_plane_id_test)
		col6 = fits.Column(name='TRK_FLAG', format='I', array=self.Glob_Si_id_test)
		col7 = fits.Column(name='STRIP_ID', format='J', array=self.Glob_Strip_id_test)
		col8 = fits.Column(name='POS', format='F20.5', array=self.Glob_pos_test)
		col9 = fits.Column(name='ZPOS', format='F20.5', array=self.Glob_zpos_test)
		col10 = fits.Column(name='E_DEP', format='F20.5', array=self.Glob_energy_dep_test)
		col11 = fits.Column(name='PAIR_FLAG', format='I', array=self.Glob_pair_flag_test)
		
		cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Creator          = Giovanni Giannella & Simone Guidotti', ext=1)

		fits.setval(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='THELSIM release  = eASTROGAM '+astrogam_version, ext=1)
		
		fits.setval(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)
		
		fits.setval(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_trig ='+str(self.N_trig), ext=1)

		fits.setval(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='ENERGY ='+self.ene_type, ext=1)

		fits.setval(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)
		
		fits.setval(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)
		
		fits.setval(self.outdir+'/L0.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = keV', ext=1)
	

    ############ L05.fits  #####################
    # IT WRITES THE FITS FILE WITH CLUSTER DATA OF THE TRACKER #

    def writing_L05(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):

		col1 = fits.Column(name='EVT_ID', format='J', array=self.Glob_event_id_cluster)	
		col2 = fits.Column(name='TRAY_ID', format='I', array=self.Glob_tray_id_cluster)
		col3 = fits.Column(name='PLANE_ID', format='I', array=self.Glob_plane_id_cluster)
		col4 = fits.Column(name='TRK_FLAG', format='I', array=self.Glob_Si_id_cluster)
		col5 = fits.Column(name='POS', format='F.20.5', array=self.Glob_pos_cluster)
		col6 = fits.Column(name='ZPOS', format='F20.5', array=self.Glob_zpos_cluster)
		col7 = fits.Column(name='E_DEP', format='F20.5', array=self.Glob_energy_dep_cluster)
		col8 = fits.Column(name='PAIR_FLAG', format='I', array=self.Glob_pair_flag_cluster)
							
		cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
			
			
		if os.path.exists(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			
		fits.setval(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Creator          = Giovanni Giannella & Simone Guidotti', ext=1)
			
		fits.setval(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='THELSIM release  = eASTROGAM '+astrogam_version, ext=1)
			
		fits.setval(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in ='+str(N_in)+'   /Number of simulated particles', ext=1)
			
		fits.setval(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_trig ='+str(self.N_trig)+'   /Number of triggering events', ext=1)
			
		fits.setval(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='ENERGY ='+self.ene_type+'   /Simulated input energy', ext=1)
			
		fits.setval(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta ='+str(theta_type)+'   /Simulated input theta angle', ext=1)
			
		fits.setval(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi ='+str(phi_type)+'   /Simulated input phi angle', ext=1)
			
		fits.setval(self.outdir+'/L0.5.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = keV', ext=1)
			

    ###### cal_raw.fits  ########
    # IT WRITES THE FITS FILE WITH RAW DATA OF THE CALORIMETER #

    def writing_cal_raw(self, N_in, part_type, theta_type, phi_type, ifile, astrogam_version):

		
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_cal)	
		col2 = fits.Column(name='VOL_ID', format='I', array=self.vol_id_cal)
		col3 = fits.Column(name='MOTH_ID', format='J', array=self.moth_id_cal)
		col4 = fits.Column(name='E_DEP', format='F20.5', array=self.energy_dep_cal)
		col5 = fits.Column(name='X_ENT', format='F20.5', array=self.x_en_cal)
		col6 = fits.Column(name='Y_ENT', format='F20.5', array=self.y_en_cal)
		col7 = fits.Column(name='Z_ENT', format='F20.5', array=self.z_en_cal)
		col8 = fits.Column(name='X_EXIT', format='F20.5', array=self.x_ex_cal)
		col9 = fits.Column(name='Y_EXIT', format='F20.5', array=self.y_ex_cal)
		col10 = fits.Column(name='Z_EXIT', format='F20.5', array=self.z_ex_cal)
		col11 = fits.Column(name='PART_ID', format='I', array=self.part_id_cal)
		col12 = fits.Column(name='TRK_ID', format='I', array=self.trk_id_cal)
		col13 = fits.Column(name='CHILD_ID', format='I', array=self.child_id_cal)
		col14 = fits.Column(name='PROC_ID', format='I', array=self.proc_id_cal)

		
		cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy ='+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Position unit = cm', ext=1)

		fits.setval(self.outdir+'/G4.RAW.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = keV', ext=1)


			
    #############  G4_Cal.fits  #############
    # IT WRITES THE FITS FILE WITH BARS DATA OF THE CALORIMETER, CONSIDERING THE FLAGGED COMPTON/PAIR EVENTS #

    def writing_G4cal(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):

		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_cal)	
		col2 = fits.Column(name='BAR_ID', format='I', array=self.bar_id_tot)
		col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=self.bar_ene_tot)
		col4 = fits.Column(name='PAIR_FLAG', format='I', array=self.pair_flag_tot_cal)
		
		cols = fits.ColDefs([col1,col2,col3,col4])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)


    ####### COMPTON/PAIR/RAYLEIGH CAL EVENTS FITS #########
    # IT WRITES THE FITS FILE WITH BARS DATA OF THE CALORIMETER, CONSIDERING THE SEPARATED COMPTON EVENTS BY PAIR PRODUCTION ONES #

    def writing_G4_cal_compton(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):
		
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_cal_compton)	
		col2 = fits.Column(name='BAR_ID', format='I', array=self.bar_id_tot_compton)
		col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=self.bar_ene_tot_compton)
		col4 = fits.Column(name='PAIR_FLAG', format='I', array=self.pair_flag_tot_cal_compton)
		
		cols = fits.ColDefs([col1,col2,col3,col4])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.CAL.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)
		

    def writing_G4_cal_pair(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):
	
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_cal_pair)	
		col2 = fits.Column(name='BAR_ID', format='I', array=self.bar_id_tot_pair)
		col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=self.bar_ene_tot_pair)
		col4 = fits.Column(name='PAIR_FLAG', format='I', array=self.pair_flag_tot_cal_pair)
		
		cols = fits.ColDefs([col1,col2,col3,col4])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.CAL.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)

    def writing_G4_cal_ray(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):
		
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_cal_ray)	
		col2 = fits.Column(name='BAR_ID', format='I', array=self.bar_id_tot_ray)
		col3 = fits.Column(name='BAR_ENERGY', format='F20.15', array=self.bar_ene_tot_ray)
		col4 = fits.Column(name='PAIR_FLAG', format='I', array=self.pair_flag_tot_cal_ray)
		
		cols = fits.ColDefs([col1,col2,col3,col4])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.CAL.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)

			
    ######### Cal Sum.fits  ##################
    # IT WRITES THE FITS FILE WITH BARS DATA OF THE CALORIMETER, CONSIDERING THE SUMMING ENERGY FOR EVERY EVENT #

    def writing_cal_sum(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):

		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_cal_sum)	
		col2 = fits.Column(name='BAR_ENERGY', format='F20.15', array=self.bar_ene_tot_sum)

		cols = fits.ColDefs([col1,col2])
		tbhdu = fits.BinTableHDU.from_columns(cols)


		if os.path.exists(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/SUM.CAL.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)


    ####### S1 data format ##############
    # ID Type Edep VolumeID X Y Z
	# - ID = Event ID
	# - Type = event flag
	# - Edep = energy deposited in the strip/calorimeter bar
	# - volume ID = 
	# - X, Y, Z = position of the center of the volume

    def writing_S1_cal(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_CAL.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_CAL.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_CAL.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_CAL.dat', 'w')		

		for j in xrange(len(self.event_id_tot_cal)):
			data.write('{:d}\t'.format(self.event_id_tot_cal[j]))
			data.write('{:d}\t'.format(self.pair_flag_tot_cal[j]))
			data.write('{:f}\t'.format(self.bar_ene_tot[j]))
			cal_unique_id = self.bar_id_tot[j] + self.cal_vol_start
			data.write('{:d}\t'.format(cal_unique_id))
			data.write('{:f}\t'.format(0.0))
			data.write('{:f}\t'.format(0.0))
			data.write('{:f}\n'.format(0.0))
			

		data.close()


    ###### AC_raw.fits  ########
    # IT WRITES THE FITS FILE WITH RAW DATA OF THE ANTICOINCIDENCE #

    def writing_ac_raw(self, N_in, part_type, theta_type, phi_type, ifile, astrogam_version):
		
		
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_ac)	
		col2 = fits.Column(name='VOL_ID', format='I', array=self.vol_id_ac)
		col3 = fits.Column(name='MOTH_ID', format='J', array=self.moth_id_ac)
		col4 = fits.Column(name='E_DEP', format='F20.5', array=self.energy_dep_ac)
		col5 = fits.Column(name='X_ENT', format='F20.5', array=self.x_en_ac)
		col6 = fits.Column(name='Y_ENT', format='F20.5', array=self.y_en_ac)
		col7 = fits.Column(name='Z_ENT', format='F20.5', array=self.z_en_ac)
		col8 = fits.Column(name='X_EXIT', format='F20.5', array=self.x_ex_ac)
		col9 = fits.Column(name='Y_EXIT', format='F20.5', array=self.y_ex_ac)
		col10 = fits.Column(name='Z_EXIT', format='F20.5', array=self.z_ex_ac)
		col11 = fits.Column(name='PART_ID', format='I', array=self.part_id_ac)
		col12 = fits.Column(name='TRK_ID', format='I', array=self.trk_id_ac)
		col13 = fits.Column(name='CHILD_ID', format='I', array=self.child_id_ac)
		col14 = fits.Column(name='PROC_ID', format='I', array=self.proc_id_ac)

		
		cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in ='+str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy ='+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta ='+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi ='+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Position unit = cm', ext=1)

		fits.setval(self.outdir+'/G4.RAW.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = keV', ext=1)




    ########## G4_AC.fits #################
    # IT WRITES THE FITS FILE WITH PANELS DATA OF THE ANTICOINCIDENCE, CONSIDERING THE FLAGGED COMPTON/PAIR/RAYLEIGH EVENTS #

    def writing_G4ac(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):

		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_ac)	
		col2 = fits.Column(name='AC_PANEL', format='A', array=self.AC_panel)
		col3 = fits.Column(name='AC_SUBPANEL', format='I', array=self.AC_subpanel)
		col4 = fits.Column(name='E_DEP', format='F20.15', array=self.energy_dep_tot_ac)
		col5 = fits.Column(name='PAIR_FLAG', format='I', array=self.pair_flag_tot_ac)

		cols = fits.ColDefs([col1,col2,col3,col4,col5])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')


		fits.setval(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.AC.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)



    ####### COMPTON/PAIR AC EVENTS FITS #########
    # IT WRITES THE FITS FILE WITH PANELS DATA OF THE ANTICOINCIDENCE, CONSIDERING THE SEPARATED COMPTON EVENTS BY PAIR PRODUCTION ONES #

    def writing_G4_ac_compton(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):
		
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_ac_compton)	
		col2 = fits.Column(name='AC_PANEL', format='A', array=self.AC_panel_compton)
		col3 = fits.Column(name='AC_SUBPANEL', format='I', array=self.AC_subpanel_compton)
		col4 = fits.Column(name='E_DEP', format='F20.15', array=self.energy_dep_tot_ac_compton)
		col5 = fits.Column(name='PAIR_FLAG', format='I', array=self.pair_flag_tot_ac_compton)

		cols = fits.ColDefs([col1,col2,col3,col4,col5])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.AC.COMPTON.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)
		

    def writing_G4_ac_pair(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):
	
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_ac_pair)	
		col2 = fits.Column(name='AC_PANEL', format='A', array=self.AC_panel_pair)
		col3 = fits.Column(name='AC_SUBPANEL', format='I', array=self.AC_subpanel_pair)
		col4 = fits.Column(name='E_DEP', format='F20.15', array=self.energy_dep_tot_ac_pair)
		col5 = fits.Column(name='PAIR_FLAG', format='I', array=self.pair_flag_tot_ac_pair)

		cols = fits.ColDefs([col1,col2,col3,col4,col5])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.AC.PAIR.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)


    def writing_G4_ac_ray(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):
		
		col1 = fits.Column(name='EVT_ID', format='I', array=self.event_id_tot_ac_ray)	
		col2 = fits.Column(name='AC_PANEL', format='A', array=self.AC_panel_ray)
		col3 = fits.Column(name='AC_SUBPANEL', format='I', array=self.AC_subpanel_ray)
		col4 = fits.Column(name='E_DEP', format='F20.15', array=self.energy_dep_tot_ac_ray)
		col5 = fits.Column(name='PAIR_FLAG', format='I', array=self.pair_flag_tot_ac_ray)

		cols = fits.ColDefs([col1,col2,col3,col4,col5])
		tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
		if os.path.exists(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits'):
			os.remove(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
			tbhdu.writeto(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')
		else:
			tbhdu.writeto(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits')

		fits.setval(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='eASTROGAM '+astrogam_version+' Geant4 simulation', ext=1)

		fits.setval(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='N_in =  ' +str(N_in), ext=1)

		fits.setval(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy     = '+self.ene_type, ext=1)

		fits.setval(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Theta     = '+str(theta_type), ext=1)

		fits.setval(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Phi     = '+str(phi_type), ext=1)

		fits.setval(self.outdir+'/G4.AC.RAYLEIGH.eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.fits', 'COMMENT', value='Energy unit = GeV', ext=1)


    ####### S1 data format ##############

    # ID Type Edep VolumeID X Y Z
    # - ID = Event ID
    # - Type = event flag
    # - Edep = energy deposited in the strip/calorimeter bar
    # - volume ID = unique ID for the volume
    # - X, Y, Z = position of the center of the volume

    def writing_S1_ac(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_AC.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_AC.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_AC.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_AC.dat', 'w')		

		for j in xrange(len(self.event_id_tot_ac)):
			data.write('{:d}\t'.format(self.event_id_tot_ac[j]))
			data.write('{:d}\t'.format(self.pair_flag_tot_ac[j]))
			data.write('{:f}\t'.format(self.energy_dep_tot_ac[j]))
			data.write('{:d}\t'.format(self.vol_id_tot_ac[j]))
			data.write('{:f}\t'.format(0.0))
			data.write('{:f}\t'.format(0.0))
			data.write('{:f}\n'.format(0.0))
			
		data.close()

	####### AA Strip ##############

	# IT WRITES THE STRIP DATA OF THE TRACKER FOR THE KALMAN INPUT #

    def writing_AA_strip(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')


		j=0
		while j < len(self.Glob_event_id_test):
			
			where_event_eq = np.where(self.Glob_event_id_test == self.Glob_event_id_test[j])
			where_event_eq = where_event_eq[0]
				
			Glob_Si_id_test_temp = self.Glob_Si_id_test[where_event_eq]
			Glob_tray_id_test_temp  = self.Glob_tray_id_test[where_event_eq]
			Glob_plane_id_test_temp  = self.Glob_plane_id_test[where_event_eq]
			Glob_Strip_id_test_temp = self.Glob_Strip_id_test[where_event_eq]
			Glob_pos_test_temp = self.Glob_pos_test[where_event_eq]
			Glob_zpos_test_temp = self.Glob_zpos_test[where_event_eq]
			Glob_energy_dep_test_temp = self.Glob_energy_dep_test[where_event_eq]
			Glob_pair_flag_test_temp = self.Glob_pair_flag_test[where_event_eq]

			# ------------------------------------
				
			# X VIEW

			r = 0

			where_x = np.where(Glob_Si_id_test_temp == 0)
			where_x = where_x[0]				
				
			if len(where_x) != 0:				
				while r < len(where_x):
					data.write('{:d}\t'.format(self.Glob_event_id_test[j]))
					data.write('{:d}\t'.format(theta_type))
					data.write('{:d}\t'.format(phi_type))
					data.write('{:s}\t'.format(self.ene_type))
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
					data.write('{:d}\t'.format(self.Glob_event_id_test[j]))
					data.write('{:d}\t'.format(theta_type))
					data.write('{:d}\t'.format(phi_type))
					data.write('{:s}\t'.format(self.ene_type))
					data.write('{:d}\t'.format(Glob_plane_id_test_temp[where_y[r]]))
					data.write('{:f}\t'.format(Glob_zpos_test_temp[where_y[r]]))
					data.write('{:d}\t'.format(1))
					data.write('{:d}\t'.format(Glob_Strip_id_test_temp[where_y[r]]))
					data.write('{:f}\t'.format(Glob_pos_test_temp[where_y[r]]))
					data.write('{:f}\t'.format(Glob_energy_dep_test_temp[where_y[r]]))
					data.write('{:d}\n'.format(Glob_pair_flag_test_temp[where_y[r]]))

					r = r + 1


			j_max = max(where_event_eq)
			j = j_max + 1


		data.close()

	####### S1 data format ##############

	# ID Type Edep VolumeID X Y Z
	# - ID = Event ID
	# - Type = event flag
	# - Edep = energy deposited in the strip/calorimeter bar
	# - volume ID = 
	# - X, Y, Z = position of the center of the volume

    def writing_S1_trk(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_TRACKER.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_TRACKER.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_TRACKER.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_S1_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'_TRACKER.dat', 'w')
			
		j=0
		while j < len(self.Glob_event_id_test):
			
			where_event_eq = np.where(self.Glob_event_id_test == self.Glob_event_id_test[j])
			where_event_eq = where_event_eq[0]		
	
			Glob_Si_id_test_temp = self.Glob_Si_id_test[where_event_eq]
			Glob_tray_id_test_temp  = self.Glob_tray_id_test[where_event_eq]
			Glob_plane_id_test_temp  = self.Glob_plane_id_test[where_event_eq]
			Glob_Strip_id_test_temp = self.Glob_Strip_id_test[where_event_eq]
			Glob_pos_test_temp = self.Glob_pos_test[where_event_eq]
			Glob_zpos_test_temp = self.Glob_zpos_test[where_event_eq]
			Glob_energy_dep_test_temp = self.Glob_energy_dep_test[where_event_eq]
			Glob_pair_flag_test_temp = self.Glob_pair_flag_test[where_event_eq]

			# ------------------------------------
				
			# X VIEW

			r = 0

			where_x = np.where(Glob_Si_id_test_temp == 0)
			where_x = where_x[0]
			
			if len(where_x) != 0:				
				while r < len(where_x):
					
					x_unique_id = Glob_Strip_id_test_temp[where_x[r]] + (self.tracker_bottom_vol_start*Glob_tray_id_test_temp[where_x[r]] + self.x_layer_id)

					data.write('{:d}\t'.format(self.Glob_event_id_test[j]))	
					data.write('{:d}\t'.format(Glob_pair_flag_test_temp[where_x[r]]))
					data.write('{:f}\t'.format(Glob_energy_dep_test_temp[where_x[r]]))
					data.write('{:d}\t'.format(x_unique_id))
					data.write('{:f}\t'.format(Glob_pos_test_temp[where_x[r]]))
					data.write('{:f}\t'.format(0.0))					
					data.write('{:f}\n'.format(Glob_zpos_test_temp[where_x[r]]))

										
					r = r + 1
					
			# ------------------------------------

			# Y VIEW

			r = 0

			where_y = np.where(Glob_Si_id_test_temp == 1)
			where_y = where_y[0]				
				
			if len(where_y) != 0:				
				while r < len(where_y):
				
					y_unique_id = Glob_Strip_id_test_temp[where_y[r]] + (self.tracker_bottom_vol_start*Glob_tray_id_test_temp[where_y[r]] + self.y_layer_id)
					
					data.write('{:d}\t'.format(self.Glob_event_id_test[j]))	
					data.write('{:d}\t'.format(Glob_pair_flag_test_temp[where_y[r]]))
					data.write('{:f}\t'.format(Glob_energy_dep_test_temp[where_y[r]]))
					data.write('{:d}\t'.format(y_unique_id))
					data.write('{:f}\t'.format(0.0))	
					data.write('{:f}\t'.format(Glob_pos_test_temp[where_y[r]]))
					data.write('{:f}\n'.format(Glob_zpos_test_temp[where_y[r]]))
					
					r = r + 1					
			
			# ------------------------------------

			j_max = max(where_event_eq)
			j = j_max + 1


		data.close()


	##### AA cluster #####

	# IT WRITES THE CLUSTER DATA OF THE TRACKER FOR THE KALMAN INPUT #

    def writing_AA_cluster(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		
		j=0
		while j < len(self.Glob_event_id_cluster):
			
			where_event_eq = np.where(self.Glob_event_id_cluster == self.Glob_event_id_cluster[j])
			where_event_eq = where_event_eq[0] 

			Glob_Si_id_cluster_temp = self.Glob_Si_id_cluster[where_event_eq]
			Glob_tray_id_cluster_temp  = self.Glob_tray_id_cluster[where_event_eq]
			Glob_plane_id_cluster_temp  = self.Glob_plane_id_cluster[where_event_eq]
			Glob_energy_dep_cluster_temp = self.Glob_energy_dep_cluster[where_event_eq]
			Glob_pos_cluster_temp = self.Glob_pos_cluster[where_event_eq]
			Glob_zpos_cluster_temp = self.Glob_zpos_cluster[where_event_eq]
			Glob_pair_flag_cluster_temp = self.Glob_pair_flag_cluster[where_event_eq]
			Glob_Strip_number_cluster_temp = self.Glob_Strip_number_cluster[where_event_eq]


			# ------------------------------------
				
			# X VIEW

			r = 0

			where_x = np.where(Glob_Si_id_cluster_temp == 0)
			where_x = where_x[0]				
				
			if len(where_x) != 0:				
				while r < len(where_x):
					data.write('{:d}\t'.format(self.Glob_event_id_cluster[j]))
					data.write('{:d}\t'.format(theta_type))
					data.write('{:d}\t'.format(phi_type))
					data.write('{:s}\t'.format(self.ene_type))
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
					data.write('{:d}\t'.format(self.Glob_event_id_cluster[j]))
					data.write('{:d}\t'.format(theta_type))
					data.write('{:d}\t'.format(phi_type))
					data.write('{:s}\t'.format(self.ene_type))
					data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_y[r]]))
					data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_y[r]]))
					data.write('{:d}\t'.format(1))
					data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_y[r]]))
					data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_y[r]]))
					data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_y[r]]))
					data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_y[r]]))

					r = r + 1


			j_max = max(where_event_eq)
			j = j_max + 1

		data.close()
	

	##### AA cluster pairs ###############

	# IT WRITES THE ONLY PAIRS CLUSTER DATA OF THE TRACKER FOR THE KALMAN INPUT #

    def writing_AA_cluster_pairs(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')

	
		totalstrips_before = 0
		j=0
		while j < len(self.Glob_event_id_cluster):
			
			where_event_eq = np.where(self.Glob_event_id_cluster == self.Glob_event_id_cluster[j])
			where_event_eq = where_event_eq[0] 

			Glob_Si_id_cluster_temp = self.Glob_Si_id_cluster[where_event_eq]
			Glob_tray_id_cluster_temp  = self.Glob_tray_id_cluster[where_event_eq]
			Glob_plane_id_cluster_temp  = self.Glob_plane_id_cluster[where_event_eq]
			Glob_energy_dep_cluster_temp = self.Glob_energy_dep_cluster[where_event_eq]
			Glob_pos_cluster_temp = self.Glob_pos_cluster[where_event_eq]
			Glob_zpos_cluster_temp = self.Glob_zpos_cluster[where_event_eq]
			Glob_pair_flag_cluster_temp = self.Glob_pair_flag_cluster[where_event_eq]
			Glob_Strip_number_cluster_temp = self.Glob_Strip_number_cluster[where_event_eq]


			# ------------------------------------
				
			# X VIEW

			r = 0

			where_x = np.where(Glob_Si_id_cluster_temp == 0)
			where_x = where_x[0]				
				
			if len(where_x) != 0:				
				while r < len(where_x):
					if Glob_pair_flag_cluster_temp[where_x[r]] == 1 or Glob_pair_flag_cluster_temp[where_x[r]] == 4:
						data.write('{:d}\t'.format(self.Glob_event_id_cluster[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(self.ene_type))
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
					if Glob_pair_flag_cluster_temp[where_y[r]] == 1 or Glob_pair_flag_cluster_temp[where_y[r]] == 4:
						data.write('{:d}\t'.format(self.Glob_event_id_cluster[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(self.ene_type))
						data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_y[r]]))
						data.write('{:d}\t'.format(1))
						data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_y[r]]))
						data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_y[r]]))
						data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_y[r]]))

					r = r + 1


			j_max = max(where_event_eq)
			j = j_max + 1

		data.close()



	######## AA cluster compton ##############

	# IT WRITES THE ONLY COMPTON CLUSTER DATA OF THE TRACKER FOR THE KALMAN INPUT #

    def writing_AA_cluster_compton(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')


		totalstrips_before = 0
		j=0
		while j < len(self.Glob_event_id_cluster):
			
			where_event_eq = np.where(self.Glob_event_id_cluster == self.Glob_event_id_cluster[j])
			where_event_eq = where_event_eq[0] 
										   
			Glob_Si_id_cluster_temp = self.Glob_Si_id_cluster[where_event_eq]
			Glob_tray_id_cluster_temp  = self.Glob_tray_id_cluster[where_event_eq]
			Glob_plane_id_cluster_temp  = self.Glob_plane_id_cluster[where_event_eq]
			Glob_energy_dep_cluster_temp = self.Glob_energy_dep_cluster[where_event_eq]
			Glob_pos_cluster_temp = self.Glob_pos_cluster[where_event_eq]
			Glob_zpos_cluster_temp = self.Glob_zpos_cluster[where_event_eq]
			Glob_pair_flag_cluster_temp = self.Glob_pair_flag_cluster[where_event_eq]
			Glob_Strip_number_cluster_temp = self.Glob_Strip_number_cluster[where_event_eq]
										   
										   
			# ------------------------------------
										   
			# X VIEW
										   
			r = 0

			where_x = np.where(Glob_Si_id_cluster_temp == 0)
			where_x = where_x[0]				
										   
			if len(where_x) != 0:
				while r < len(where_x):
					if Glob_pair_flag_cluster_temp[where_x[r]] == 2 or Glob_pair_flag_cluster_temp[where_x[r]] == 5:
						data.write('{:d}\t'.format(self.Glob_event_id_cluster[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(self.ene_type))
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
					if Glob_pair_flag_cluster_temp[where_y[r]] == 2 or Glob_pair_flag_cluster_temp[where_y[r]] == 5:
						data.write('{:d}\t'.format(self.Glob_event_id_cluster[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(self.ene_type))
						data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_y[r]]))
						data.write('{:d}\t'.format(1))
						data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_y[r]]))
						data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_y[r]]))
						data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_y[r]]))
															 
					r = r + 1
															 
															 
			j_max = max(where_event_eq)
			j = j_max + 1

		data.close()
	
	######## AA cluster rayleigh ##############

	# IT WRITES THE ONLY RAYLEIGH CLUSTER DATA OF THE TRACKER FOR THE KALMAN INPUT #

    def writing_AA_cluster_rayleigh(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_CLUSTER_RAYLEIGH_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')


		totalstrips_before = 0
		j=0
		while j < len(self.Glob_event_id_cluster):
			
			where_event_eq = np.where(self.Glob_event_id_cluster == self.Glob_event_id_cluster[j])
			where_event_eq = where_event_eq[0] 
										   
			Glob_Si_id_cluster_temp = self.Glob_Si_id_cluster[where_event_eq]
			Glob_tray_id_cluster_temp  = self.Glob_tray_id_cluster[where_event_eq]
			Glob_plane_id_cluster_temp  = self.Glob_plane_id_cluster[where_event_eq]
			Glob_energy_dep_cluster_temp = self.Glob_energy_dep_cluster[where_event_eq]
			Glob_pos_cluster_temp = self.Glob_pos_cluster[where_event_eq]
			Glob_zpos_cluster_temp = self.Glob_zpos_cluster[where_event_eq]
			Glob_pair_flag_cluster_temp = self.Glob_pair_flag_cluster[where_event_eq]
			Glob_Strip_number_cluster_temp = self.Glob_Strip_number_cluster[where_event_eq]
										   
										   
			# ------------------------------------
										   
			# X VIEW
										   
			r = 0

			where_x = np.where(Glob_Si_id_cluster_temp == 0)
			where_x = where_x[0]				
										   
			if len(where_x) != 0:
				while r < len(where_x):
					if Glob_pair_flag_cluster_temp[where_x[r]] == 3 or Glob_pair_flag_cluster_temp[where_x[r]] == 6:
						data.write('{:d}\t'.format(self.Glob_event_id_cluster[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(self.ene_type))
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
					if Glob_pair_flag_cluster_temp[where_y[r]] == 3 or Glob_pair_flag_cluster_temp[where_y[r]] == 6:
						data.write('{:d}\t'.format(self.Glob_event_id_cluster[j]))
						data.write('{:d}\t'.format(theta_type))
						data.write('{:d}\t'.format(phi_type))
						data.write('{:s}\t'.format(self.ene_type))
						data.write('{:d}\t'.format(Glob_plane_id_cluster_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_zpos_cluster_temp[where_y[r]]))
						data.write('{:d}\t'.format(1))
						data.write('{:f}\t'.format(Glob_pos_cluster_temp[where_y[r]]))
						data.write('{:f}\t'.format(Glob_energy_dep_cluster_temp[where_y[r]]))
						data.write('{:d}\t'.format(Glob_Strip_number_cluster_temp[where_y[r]]))
						data.write('{:d}\n'.format(Glob_pair_flag_cluster_temp[where_y[r]]))
															 
					r = r + 1
															 
															 
			j_max = max(where_event_eq)
			j = j_max + 1

		data.close()
	
	######### AA Fake #########

	# IT WRITES THE RAW DATA OF THE TRACKER WHEN THE STRIPS AREN'T ACTIVED #

    def writing_AA_fake(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):

		if os.path.exists(self.outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data = open(self.outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data = open(self.outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'.'+self.py_name+'.'+self.sim_name+'.'+self.stripname+'.'+self.sname+'.'+str(N_in)+part_type+'.'+self.ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
	
		j = 0
		while j < len(self.event_id_tr):
			
			where_event_eq = np.where(self.event_id_tr == self.event_id_tr[j])
			where_event_eq = where_event_eq[0]

			plane_id_temp = self.plane_id[where_event_eq]
			Cluster_x_temp  = self.x_pos[where_event_eq]
			Cluster_y_temp  = self.y_pos[where_event_eq]
			Cluster_z_temp  = self.z_pos[where_event_eq]
			e_dep_x_temp  = (self.en_dep_tr[where_event_eq])/2.
			e_dep_y_temp  = (self.en_dep_tr[where_event_eq])/2.
			child_temp = self.child_id_tr[where_event_eq]
			proc_temp = self.proc_id_tr[where_event_eq]

			# ------------------------------------
				
			# X VIEW

			r = 0
			while r < len(Cluster_x_temp):

				if e_dep_x_temp[r] > self.E_th:
					data.write('{:d}\t'.format(self.event_id_tr[j]))
					data.write('{:d}\t'.format(theta_type))
					data.write('{:d}\t'.format(phi_type))
					data.write('{:s}\t'.format(self.ene_type))
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
				if e_dep_y_temp[r] > self.E_th:
					data.write('{:d}\t'.format(self.event_id_tr[j]))
					data.write('{:d}\t'.format(theta_type))
					data.write('{:d}\t'.format(phi_type))
					data.write('{:s}\t'.format(self.ene_type))
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


    ############# LOADING GEOM ################

	# IT LOADS THE LUT GENERATED BY eASTROGAM_GEOMETRY

    def loading_geom(self, astrogam_version):

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


		self.Arch_vol_id_x_top = Arch_vol_id_x_top 
		self.Arch_moth_id_x_top = Arch_moth_id_x_top 
		self.Arch_Strip_id_x_top = Arch_Strip_id_x_top 
		self.Arch_Si_id_x_top = Arch_Si_id_x_top 
		self.Arch_tray_id_x_top = Arch_tray_id_x_top 
		self.Arch_plane_id_x_top = Arch_plane_id_x_top 
		self.Arch_xpos_x_top = Arch_xpos_x_top 
		self.Arch_zpos_x_top = Arch_zpos_x_top 
		self.Arch_energy_dep_x_top = Arch_energy_dep_x_top 
		self.Arch_pair_flag_x_top = Arch_pair_flag_x_top 

		self.Arch_vol_id_y_top = Arch_vol_id_y_top 
		self.Arch_moth_id_y_top = Arch_moth_id_y_top 
		self.Arch_Strip_id_y_top = Arch_Strip_id_y_top 
		self.Arch_Si_id_y_top = Arch_Si_id_y_top 
		self.Arch_tray_id_y_top = Arch_tray_id_y_top 
		self.Arch_plane_id_y_top = Arch_plane_id_y_top 
		self.Arch_ypos_y_top = Arch_ypos_y_top  
		self.Arch_zpos_y_top = Arch_zpos_y_top 
		self.Arch_energy_dep_y_top = Arch_energy_dep_y_top 
		self.Arch_pair_flag_y_top = Arch_pair_flag_y_top
	
     
######### FLAG EVENTS ##############

	# IT FLAGS THE PRIMARY EVENT, SEPARATING COMPTON FROM PAIR FROM RAYLEIGH FROM OTHER #
	# IT SUMS THE ENERGY OF THE SAME EVENT FOR EVERY VOLUME #
	# IT FLAGS THE SUMMED EVENT IF ONE HITS ARE IN THE SAME VOLUME OR IT FLAGS THE SUMMED EVENT IF THE HITS AREN'T IN THE SAME VOLUME  #
	# IT PUTS EVERYTHING IN THE SAME ARRAY #

    def flag_events(self):

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
		gtime_tot = []
			
		j = 0
		while j < len(self.event_id_tr):          
		
			where_event_eq = np.where(self.event_id_tr == self.event_id_tr[j])
			where_event_eq = where_event_eq[0]
				
			vol_id_temp = self.vol_id_tr[where_event_eq]
			moth_id_temp  = self.moth_id_tr[where_event_eq]
			Strip_id_x_temp  = self.Strip_id_x[where_event_eq]
			Strip_id_y_temp  = self.Strip_id_y[where_event_eq]
			tray_id_temp  = self.tray_id[where_event_eq]
			plane_id_temp  = self.plane_id[where_event_eq]
			energy_dep_temp = self.en_dep_tr[where_event_eq]
			child_id_temp = self.child_id_tr[where_event_eq]
			proc_id_temp = self.proc_id_tr[where_event_eq]
			gtime_temp = self.gtime_ent_tr[where_event_eq]
			trk_id_temp = self.trk_id_tr[where_event_eq]
			theta_ent_id_temp = self.theta_ent_tr[where_event_eq]
			phi_ent_id_temp = self.phi_ent_tr[where_event_eq]
			theta_exit_id_temp = self.theta_exit_tr[where_event_eq]
			phi_exit_id_temp = self.phi_exit_tr[where_event_eq]

			where_pair = np.where((child_id_temp == 1) & (proc_id_temp == 7) & (trk_id_temp <= 3))
			where_pair = where_pair[0]
			where_compton = np.where((child_id_temp == 0) & (theta_ent_id_temp != theta_exit_id_temp) & (phi_ent_id_temp != phi_exit_id_temp) & (energy_dep_temp > 0.))
			where_compton = where_compton[0]			
			#where_compton = np.where((child_id_temp == 1) & (proc_id_temp == 3) & (trk_id_temp <= 2))
			#where_compton = where_compton[0]
			where_ray = np.where((child_id_temp == 0) & (theta_ent_id_temp != theta_exit_id_temp) & (phi_ent_id_temp != phi_exit_id_temp) & (energy_dep_temp == 0.))
			where_ray = where_ray[0]
			#where_ray = np.where((child_id_temp == 1) & (proc_id_temp == 9) & (trk_id_temp == 1))
			#where_ray = where_ray[0]
			
			ispair = 0
			iscompton = 0
			isray = 0
			isother = 0
			gtime_pair = [10**9]
			gtime_compton = [10**9]
			gtime_ray = [10**9]
			
			if len(where_pair) != 0:
				gtime_pair = gtime_temp[where_pair]
			else:
				gtime_pair = [10**9]

			if len(where_compton) != 0:
				gtime_compton = gtime_temp[where_compton]
			else:
				gtime_compton = [10**9]

			if len(where_ray) != 0:
				gtime_ray = gtime_temp[where_ray]
			else:
				gtime_ray = [10**9]
			
			if ((gtime_pair[0] == 10**9) and (gtime_compton[0] == 10**9) and (gtime_ray[0] == 10**9)):
				isother = 1
			else:
				gtime_array = [gtime_pair[0], gtime_compton[0], gtime_ray[0]]
				proc_index = np.argmin(gtime_array)
				if proc_index == 0: ispair = 1
				if proc_index == 1: iscompton = 1
				if proc_index == 2: isray = 1

		
			r = 0												
			while 1:
						
				ispair_vol = 0
				iscompton_vol = 0
				isray_vol = 0
				isother_vol = 0
				isprimary_vol = 0
	
				where_vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]))										
				where_vol_eq = where_vol_eq[0]	

				where_other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]))
				where_other_vol = where_other_vol[0]
		
				# summing the energy of the same event 
				e_dep_temp_old = np.sum(energy_dep_temp[where_vol_eq])
				event_id_tot_old = self.event_id_tr[j]
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
					
				#Searching for Pair/Compton events of secondary particles generated by the primary
				#if another process is involved, the event is flagged as 0
				#if one of hits in the same volume is a pair the summed event is flagged as 1
				#if one of hits in the same volume is a compton the summed event is flagged as 2
				#if one of hits in the same volume is a rayleigh the summed event is flagged as 3


				all_child = child_id_temp[where_vol_eq]
				all_trk = trk_id_temp[where_vol_eq]
				all_proc = proc_id_temp[where_vol_eq]
				all_gtime = gtime_temp[where_vol_eq]
				all_theta_ent = theta_ent_id_temp[where_vol_eq]
				all_phi_ent = phi_ent_id_temp[where_vol_eq]
				all_theta_exit = theta_exit_id_temp[where_vol_eq]
				all_phi_exit = phi_exit_id_temp[where_vol_eq]
				

				where_pair_vol = np.where((all_child == 1) & (all_proc == 7) & (all_trk <= 3))
				where_pair_vol = where_pair_vol[0]

				where_compton_vol = np.where((all_child == 1) & (all_proc == 3) & (all_trk <= 2))
				where_compton_vol = where_compton_vol[0]

				#where_ray_vol = np.where((all_child == 1) & (all_proc == 9) & (all_trk == 1))
				#where_ray_vol = where_ray_vol[0]


				if len(where_pair_vol) != 0 and ispair == 1:
					pair_flag_tot_old = 1
					ispair_vol = 1

					pair_flag_tot.append(pair_flag_tot_old)	
	 

				if len(where_compton_vol) != 0 and iscompton == 1:
					pair_flag_tot_old = 2
					iscompton_vol = 1
				
					pair_flag_tot.append(pair_flag_tot_old)	

				#if len(where_ray_vol) != 0 and isray == 1:
				if ((len(where_pair_vol) == 0) and (len(where_compton_vol) == 0)):
					if isray == 1:
						pair_flag_tot_old = 3
						isray_vol = 1
						
						pair_flag_tot.append(pair_flag_tot_old)	
					else:				
						pair_flag_tot_old = 0
						pair_flag_tot.append(pair_flag_tot_old)
				else:
					if ((len(where_pair_vol) != 0 and ispair == 0) | (len(where_compton_vol) != 0 and iscompton == 0)):
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
					trk_id_temp = trk_id_temp[where_other_vol]
					child_id_temp = child_id_temp[where_other_vol]
					proc_id_temp = proc_id_temp[where_other_vol]
					gtime_temp = gtime_temp[where_other_vol]
					
				else:
					break
			
			j_max = max(where_event_eq)
			j = j_max + 1
			

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
		
		self.event_id_tot_tr_raw = event_id_tot
		self.pair_flag_tot_tr_raw = pair_flag_tot

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

		self.event_id_tot_temp = event_id_tot_temp 
		self.vol_id_tot_temp = vol_id_tot_temp 
		self.moth_id_tot_temp = moth_id_tot_temp 
		self.Strip_id_tot_temp = Strip_id_tot_temp 
		self.Si_id_tot_temp = Si_id_tot_temp 
		self.tray_id_tot_temp = tray_id_tot_temp 
		self.plane_id_tot_temp = plane_id_tot_temp 
		self.energy_dep_tot_temp = energy_dep_tot_temp 
		self.pair_flag_tot_temp = pair_flag_tot_temp


    def acap(self, astrogam_version):
		event_id_tot_list = []
		vol_id_tot_list = []
		moth_id_tot_list = []
		Strip_id_tot_list = []
		Si_id_tot_list = []
		tray_id_tot_list = []
		plane_id_tot_list = []
		#energy_dep_tot_list = []
		pair_flag_tot_list = []
		
		energy = np.zeros(len(self.energy_dep_tot_temp))

		j = 0
		while j < len(self.energy_dep_tot_temp):	
			event_id_tot = self.event_id_tot_temp[j]
			vol_id_tot = self.vol_id_tot_temp[j]
			moth_id_tot = self.moth_id_tot_temp[j]
			Strip_id_tot = self.Strip_id_tot_temp[j]
			Si_id_tot = self.Si_id_tot_temp[j]
			tray_id_tot = self.tray_id_tot_temp[j]
			plane_id_tot = self.plane_id_tot_temp[j]
			pair_flag_tot = self.pair_flag_tot_temp[j]
			
			event_id_tot_list.append(event_id_tot)
			vol_id_tot_list.append(vol_id_tot)
			moth_id_tot_list.append(moth_id_tot)
			Strip_id_tot_list.append(Strip_id_tot)
			Si_id_tot_list.append(Si_id_tot)
			tray_id_tot_list.append(tray_id_tot)
			plane_id_tot_list.append(plane_id_tot)
			pair_flag_tot_list.append(pair_flag_tot)
			
			if astrogam_version == 'V1.0' or astrogam_version == 'V1.1' or astrogam_version == 'V2.0' or astrogam_version=='V10.0':
				if (j == 0 and self.energy_dep_tot_temp[j] > 0.):
					energy[j]=0.68*self.energy_dep_tot_temp[j]
					energy[j+1]=energy[j+1]+0.17*self.energy_dep_tot_temp[j]
					energy[j+2]=energy[j+2]+0.03*self.energy_dep_tot_temp[j]
				if (j == 1 and self.energy_dep_tot_temp[j] > 0.):
					energy[j-1]=energy[j-1]+0.17*self.energy_dep_tot_temp[j]
					energy[j]=0.68*self.energy_dep_tot_temp[j]
					energy[j+1]=energy[j+1]+0.17*self.energy_dep_tot_temp[j]
					energy[j+2]=energy[j+2]+0.03*self.energy_dep_tot_temp[j]
				if (j > 1 and j < len(self.energy_dep_tot_temp)-2 and self.energy_dep_tot_temp[j] > 0.):
					energy[j-2]=energy[j-2]+0.03*self.energy_dep_tot_temp[j]
					energy[j-1]=energy[j-1]+0.17*self.energy_dep_tot_temp[j]
					energy[j]=0.68*self.energy_dep_tot_temp[j]
					energy[j+1]=energy[j+1]+0.17*self.energy_dep_tot_temp[j]
					energy[j+2]=energy[j+2]+0.03*self.energy_dep_tot_temp[j]
				if (j == len(self.energy_dep_tot_temp)-2 and self.energy_dep_tot_temp[j] > 0.):
					energy[j-2]=energy[j-2]+0.03*self.energy_dep_tot_temp[j]
					energy[j-1]=energy[j-1]+0.17*self.energy_dep_tot_temp[j]
					energy[j]=0.68*self.energy_dep_tot_temp[j]
					energy[j+1]=energy[j+1]+0.17*self.energy_dep_tot_temp[j]
				if (j == len(self.energy_dep_tot_temp)-1 and self.energy_dep_tot_temp[j] > 0.):
					energy[j-2]=energy[j-2]+0.03*self.energy_dep_tot_temp[j]
					energy[j-1]=energy[j-1]+0.17*self.energy_dep_tot_temp[j]
					energy[j]=0.68*self.energy_dep_tot_temp[j]
				
				else:
					energy[j] = self.energy_dep_tot_temp[j]

				if j == len(self.energy_dep_tot_temp)-1:
					break
				else:
					j = j + 1
					
			if astrogam_version == 'V1.2':
				if (j == 0 and self.energy_dep_tot_temp[j] > 0.):
					energy[j]=0.94*self.energy_dep_tot_temp[j]
					energy[j+1]=energy[j+1]+0.03*self.energy_dep_tot_temp[j]
				if (j == 1 and self.energy_dep_tot_temp[j] > 0.):
					energy[j-1]=energy[j-1]+0.03*self.energy_dep_tot_temp[j]
					energy[j]=0.94*self.energy_dep_tot_temp[j]
					energy[j+1]=energy[j+1]+0.03*self.energy_dep_tot_temp[j]
				if (j > 1 and j < len(self.energy_dep_tot_temp)-2 and self.energy_dep_tot_temp[j] > 0.):
					energy[j-1]=energy[j-1]+0.03*self.energy_dep_tot_temp[j]
					energy[j]=0.94*self.energy_dep_tot_temp[j]
					energy[j+1]=energy[j+1]+0.03*self.energy_dep_tot_temp[j]
				if (j == len(self.energy_dep_tot_temp)-2 and self.energy_dep_tot_temp[j] > 0.):
					energy[j-1]=energy[j-1]+0.03*self.energy_dep_tot_temp[j]
					energy[j]=0.94*self.energy_dep_tot_temp[j]
					energy[j+1]=energy[j+1]+0.03*self.energy_dep_tot_temp[j]
				if (j == len(self.energy_dep_tot_temp)-1 and self.energy_dep_tot_temp[j] > 0.):
					energy[j-1]=energy[j-1]+0.03*self.energy_dep_tot_temp[j]
					energy[j]=0.94*self.energy_dep_tot_temp[j]
				
				else:
					energy[j] = self.energy_dep_tot_temp[j]

				if j == len(self.energy_dep_tot_temp)-1:
					break
				else:
					j = j + 1		
		event_id_tot_temp = np.array(event_id_tot_list)
		vol_id_tot_temp = np.array(vol_id_tot_list)
		moth_id_tot_temp = np.array(moth_id_tot_list)
		Strip_id_tot_temp = np.array(Strip_id_tot_list)
		Si_id_tot_temp = np.array(Si_id_tot_list)
		tray_id_tot_temp = np.array(tray_id_tot_list)
		plane_id_tot_temp = np.array(plane_id_tot_list)
		energy_dep_tot_temp = np.array(energy)
		pair_flag_tot_temp = np.array(pair_flag_tot_list)
		
		self.event_id_tot_temp = event_id_tot_temp 
		self.vol_id_tot_temp = vol_id_tot_temp 
		self.moth_id_tot_temp = moth_id_tot_temp
		self.Strip_id_tot_temp = Strip_id_tot_temp 
		self.Si_id_tot_temp = Si_id_tot_temp 
		self.tray_id_tot_temp = tray_id_tot_temp 
		self.plane_id_tot_temp = plane_id_tot_temp
		self.energy_dep_tot_temp = energy_dep_tot_temp 
		self.pair_flag_tot_temp = pair_flag_tot_temp


    def noise(self):
		### 1keV di noise per ogni strip ---> *1000 (eV) ---> /2.35 (FWHM) ---> /3.6 (conv) ----> 118 e- 
		keV = 1.
		noise = keV*1000/2.35/3.6
			
		event_id_tot_list = []
		vol_id_tot_list = []
		moth_id_tot_list = []
		Strip_id_tot_list = []
		Si_id_tot_list = []
		tray_id_tot_list = []
		plane_id_tot_list = []
		#energy_dep_tot_list = []
		pair_flag_tot_list = []
		
		energy = np.zeros(len(self.energy_dep_tot_temp))
		
		j = 0
		while j < len(self.energy_dep_tot_temp):
			event_id_tot = self.event_id_tot_temp[j]
			vol_id_tot = self.vol_id_tot_temp[j]
			moth_id_tot = self.moth_id_tot_temp[j]
			Strip_id_tot = self.Strip_id_tot_temp[j]
			Si_id_tot = self.Si_id_tot_temp[j]
			tray_id_tot = self.tray_id_tot_temp[j]
			plane_id_tot = self.plane_id_tot_temp[j]
			pair_flag_tot = self.pair_flag_tot_temp[j]
			
			event_id_tot_list.append(event_id_tot)
			vol_id_tot_list.append(vol_id_tot)
			moth_id_tot_list.append(moth_id_tot)
			Strip_id_tot_list.append(Strip_id_tot)
			Si_id_tot_list.append(Si_id_tot)
			tray_id_tot_list.append(tray_id_tot)
			plane_id_tot_list.append(plane_id_tot)
			pair_flag_tot_list.append(pair_flag_tot)
		
			energy[j] = (((self.energy_dep_tot_temp[j]*1000.)/2.35)/3.6 + noise)*3.6*2.35/1000.
			
			j = j + 1
		
		event_id_tot_temp = np.array(event_id_tot_list)
		vol_id_tot_temp = np.array(vol_id_tot_list)
		moth_id_tot_temp = np.array(moth_id_tot_list)
		Strip_id_tot_temp = np.array(Strip_id_tot_list)
		Si_id_tot_temp = np.array(Si_id_tot_list)
		tray_id_tot_temp = np.array(tray_id_tot_list)
		plane_id_tot_temp = np.array(plane_id_tot_list)
		energy_dep_tot_temp = np.array(energy)
		pair_flag_tot_temp = np.array(pair_flag_tot_list)
		
		self.event_id_tot_temp = event_id_tot_temp 
		self.vol_id_tot_temp = vol_id_tot_temp 
		self.moth_id_tot_temp = moth_id_tot_temp 
		self.Strip_id_tot_temp = Strip_id_tot_temp 
		self.Si_id_tot_temp = Si_id_tot_temp
		self.tray_id_tot_temp = tray_id_tot_temp
		self.plane_id_tot_temp = plane_id_tot_temp
		self.energy_dep_tot_temp = energy_dep_tot_temp
		self.pair_flag_tot_temp = pair_flag_tot_temp



######### SUMMING ENERGY AND ENERGY THRESHOLD ###########

	# Summing the energy along the strip 

    def summing_energy(self, astrogam_version, N_in, part_type, theta_type, phi_type, ifile):

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
		while j < len(self.event_id_tot_temp):
		           
			where_event_eq = np.where(self.event_id_tot_temp == self.event_id_tot_temp[j])
			where_event_eq = where_event_eq[0]


			vol_id_temp = self.vol_id_tot_temp[where_event_eq]
			moth_id_temp = self.moth_id_tot_temp[where_event_eq]
			Strip_id_temp = self.Strip_id_tot_temp[where_event_eq]
			Si_id_temp = self.Si_id_tot_temp[where_event_eq]
			tray_id_temp = self.tray_id_tot_temp[where_event_eq]
			plane_id_temp = self.plane_id_tot_temp[where_event_eq]
			energy_dep_temp = self.energy_dep_tot_temp[where_event_eq]
			pair_flag_temp = self.pair_flag_tot_temp[where_event_eq]
				

				
			r = 0												
			while 1:

				where_vol_eq = np.where((vol_id_temp == vol_id_temp[r]) & (moth_id_temp == moth_id_temp[r]) & (Si_id_temp == 0))	
				where_vol_eq = where_vol_eq[0]
					
				where_other_vol = np.where((vol_id_temp != vol_id_temp[r]) | (moth_id_temp != moth_id_temp[r]) | (Si_id_temp != 0))
				where_other_vol = where_other_vol[0]

				if len(where_vol_eq) != 0:
					e_dep_temp_old = np.sum(energy_dep_temp[where_vol_eq])
					event_id_tot_old = self.event_id_tot_temp[j]
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
					event_id_tot_old = self.event_id_tot_temp[j]
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

			j_max = max(where_event_eq)
			j = j_max + 1
			
			
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

		self.event_id_tot = event_id_tot
		self.vol_id_tot = vol_id_tot
		self.moth_id_tot = moth_id_tot
		self.Strip_id_tot = Strip_id_tot
		self.Si_id_tot = Si_id_tot
		self.tray_id_tot = tray_id_tot
		self.plane_id_tot = plane_id_tot
		self.energy_dep_tot = energy_dep_tot
		self.pair_flag_tot = pair_flag_tot

		# apply the energy thresold

    def energy_threshold(self):
			
		where_eth = np.where(self.energy_dep_tot >= self.E_th)
		where_eth = where_eth[0]

		event_id_tot = self.event_id_tot[where_eth]
		vol_id_tot = self.vol_id_tot[where_eth]
		moth_id_tot = self.moth_id_tot[where_eth]
		Strip_id_tot = self.Strip_id_tot[where_eth]
		Si_id_tot = self.Si_id_tot[where_eth]
		tray_id_tot = self.tray_id_tot[where_eth]
		plane_id_tot = self.plane_id_tot[where_eth]
		energy_dep_tot = self.energy_dep_tot[where_eth]
		pair_flag_tot = self.pair_flag_tot[where_eth]
		
		self.event_id_tot = event_id_tot 
		self.vol_id_tot = vol_id_tot 
		self.moth_id_tot = moth_id_tot 
		self.Strip_id_tot = Strip_id_tot 
		self.Si_id_tot = Si_id_tot 
		self.tray_id_tot = tray_id_tot 
		self.plane_id_tot = plane_id_tot 
		self.energy_dep_tot = energy_dep_tot 
		self.pair_flag_tot = pair_flag_tot		


##### INDICI UNICI ##############

	# IT CREATES AN ARRAY IN WHICH THE EVENT ID IS WRITTEN ONLY ONE TIME #	

    def index_uniq(self):
									
		event_id_tot_index = []							
									
		for b in range(len(self.event_id_tot)):						
			event_id_tot_index_old = b						
			event_id_tot_index.append(event_id_tot_index_old)				
									
		event_id_tot_uniq_index = []						
									
		b = 1								
		while b < len(self.event_id_tot):						
									
			if self.event_id_tot[b] != self.event_id_tot[b-1]:					
				event_id_tot_uniq_index.append(event_id_tot_index[b-1])		
									
			if b == len(self.event_id_tot)-1:					
				event_id_tot_uniq_index.append(event_id_tot_index[b])			
									
			b = b + 1						
	
		N_trig = len(event_id_tot_uniq_index)

		event_array = self.event_id_tot[event_id_tot_uniq_index]

		self.event_array = event_array
		self.N_trig = N_trig							
  		
#### STRIP ANALYSIS ########

	# IT ANALYZES EVERY STRIP #
	
    def strip_analysis(self):


		#Total number of strips
		Total_vol_x_top = (self.N_tray)*self.N_strip
		Total_vol_y_top = (self.N_tray)*self.N_strip

		Glob_event_id_x_top = np.zeros((Total_vol_x_top, self.N_trig), dtype = np.int64)
		Glob_vol_id_x_top = np.zeros((Total_vol_x_top, self.N_trig), dtype = np.int64)
		Glob_moth_id_x_top = np.zeros((Total_vol_x_top, self.N_trig), dtype = np.int64)
		Glob_Strip_id_x_top = np.zeros((Total_vol_x_top, self.N_trig), dtype = np.int64)
		Glob_Si_id_x_top = np.zeros((Total_vol_x_top, self.N_trig), dtype = np.int64)
		Glob_tray_id_x_top = np.zeros((Total_vol_x_top, self.N_trig), dtype = np.int64)
		Glob_plane_id_x_top = np.zeros((Total_vol_x_top, self.N_trig), dtype = np.int64)
		Glob_xpos_x_top = np.zeros((Total_vol_x_top, self.N_trig))
		Glob_zpos_x_top = np.zeros((Total_vol_x_top, self.N_trig))
		Glob_energy_dep_x_top = np.zeros((Total_vol_x_top, self.N_trig))
		Glob_pair_flag_x_top = np.zeros((Total_vol_x_top, self.N_trig))

		Glob_event_id_y_top = np.zeros((Total_vol_y_top, self.N_trig), dtype = np.int64)
		Glob_vol_id_y_top = np.zeros((Total_vol_y_top, self.N_trig), dtype = np.int64)
		Glob_moth_id_y_top = np.zeros((Total_vol_y_top, self.N_trig), dtype = np.int64)
		Glob_Strip_id_y_top = np.zeros((Total_vol_y_top, self.N_trig), dtype = np.int64)
		Glob_Si_id_y_top = np.zeros((Total_vol_y_top, self.N_trig), dtype = np.int64)
		Glob_tray_id_y_top = np.zeros((Total_vol_y_top, self.N_trig), dtype = np.int64)
		Glob_plane_id_y_top = np.zeros((Total_vol_y_top, self.N_trig), dtype = np.int64)
		Glob_ypos_y_top = np.zeros((Total_vol_y_top, self.N_trig))
		Glob_zpos_y_top = np.zeros((Total_vol_y_top, self.N_trig))
		Glob_energy_dep_y_top = np.zeros((Total_vol_y_top, self.N_trig))
		Glob_pair_flag_y_top = np.zeros((Total_vol_y_top, self.N_trig))



		for i in range(self.N_trig):
			Glob_vol_id_x_top[:,i] = self.Arch_vol_id_x_top
			Glob_moth_id_x_top[:,i] = self.Arch_moth_id_x_top
			Glob_Strip_id_x_top[:,i] = self.Arch_Strip_id_x_top
			Glob_Si_id_x_top[:,i] = self.Arch_Si_id_x_top
			Glob_tray_id_x_top[:,i] = self.Arch_tray_id_x_top
			Glob_plane_id_x_top[:,i] = self.Arch_plane_id_x_top
			Glob_xpos_x_top[:,i] = self.Arch_xpos_x_top
			Glob_zpos_x_top[:,i] = self.Arch_zpos_x_top
			Glob_energy_dep_x_top[:,i] = self.Arch_energy_dep_x_top
			Glob_pair_flag_x_top[:,i] = self.Arch_pair_flag_x_top

			Glob_vol_id_y_top[:,i] = self.Arch_vol_id_y_top
			Glob_moth_id_y_top[:,i] = self.Arch_moth_id_y_top
			Glob_Strip_id_y_top[:,i] = self.Arch_Strip_id_y_top
			Glob_Si_id_y_top[:,i] = self.Arch_Si_id_y_top
			Glob_tray_id_y_top[:,i] = self.Arch_tray_id_y_top
			Glob_plane_id_y_top[:,i] = self.Arch_plane_id_y_top
			Glob_ypos_y_top[:,i] = self.Arch_ypos_y_top
			Glob_zpos_y_top[:,i] = self.Arch_zpos_y_top
			Glob_energy_dep_y_top[:,i] = self.Arch_energy_dep_y_top
			Glob_pair_flag_y_top[:,i] = self.Arch_pair_flag_y_top


		j = 0
		N_ev = 0
		while j < len(self.event_id_tot):
					
			where_event_eq = np.where(self.event_id_tot == self.event_id_tot[j])
			where_event_eq = where_event_eq[0]
				
			event_id_temp = self.event_id_tot[where_event_eq]
			vol_id_temp = self.vol_id_tot[where_event_eq]
			moth_id_temp = self.moth_id_tot[where_event_eq]
			Strip_id_temp = self.Strip_id_tot[where_event_eq]
			Si_id_temp = self.Si_id_tot[where_event_eq]
			tray_id_temp = self.tray_id_tot[where_event_eq]
			plane_id_temp = self.plane_id_tot[where_event_eq]
			energy_dep_temp = self.energy_dep_tot[where_event_eq]
			pair_flag_temp = self.pair_flag_tot[where_event_eq]			


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

					

			j_max = max(where_event_eq)			
			j = j_max + 1
			
			if j < len(self.event_id_tot):
				N_ev = N_ev + 1
		
		self.Glob_vol_id_x_top = Glob_vol_id_x_top 
		self.Glob_moth_id_x_top = Glob_moth_id_x_top 
		self.Glob_Strip_id_x_top = Glob_Strip_id_x_top 
		self.Glob_Si_id_x_top = Glob_Si_id_x_top 
		self.Glob_tray_id_x_top = Glob_tray_id_x_top 
		self.Glob_plane_id_x_top = Glob_plane_id_x_top 
		self.Glob_xpos_x_top = Glob_xpos_x_top 
		self.Glob_zpos_x_top = Glob_zpos_x_top 
		self.Glob_energy_dep_x_top = Glob_energy_dep_x_top 
		self.Glob_pair_flag_x_top = Glob_pair_flag_x_top 
		
		self.Glob_vol_id_y_top = Glob_vol_id_y_top 
		self.Glob_moth_id_y_top = Glob_moth_id_y_top 
		self.Glob_Strip_id_y_top = Glob_Strip_id_y_top 
		self.Glob_Si_id_y_top = Glob_Si_id_y_top 
		self.Glob_tray_id_y_top = Glob_tray_id_y_top 
		self.Glob_plane_id_y_top = Glob_plane_id_y_top 
		self.Glob_ypos_y_top = Glob_ypos_y_top 
		self.Glob_zpos_y_top = Glob_zpos_y_top
		self.Glob_energy_dep_y_top = Glob_energy_dep_y_top 
		self.Glob_pair_flag_y_top = Glob_pair_flag_y_top 
		self.N_ev = N_ev	

		
######## BUILD L0 ##########

	# IT CONSIDERS THE STRIP WHERE THE ENERGY IS DEPOSITED TO BUILD THE L0 OUTPUT #

    def build_L0(self):

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

		for j in range(self.N_trig):

			where_test_x = np.where(self.Glob_energy_dep_x_top[:,j] > 0.)
			where_test_x = where_test_x[0]
				
			if len(where_test_x) != 0:
				Glob_vol_id_x_test_temp = self.Glob_vol_id_x_top[where_test_x,j]
				Glob_moth_id_x_test_temp = self.Glob_moth_id_x_top[where_test_x,j]
				Glob_Strip_id_x_test_temp = self.Glob_Strip_id_x_top[where_test_x,j]
				Glob_Si_id_x_test_temp = self.Glob_Si_id_x_top[where_test_x,j]
				Glob_tray_id_x_test_temp = self.Glob_tray_id_x_top[where_test_x,j]
				Glob_plane_id_x_test_temp = self.Glob_plane_id_x_top[where_test_x,j]
				Glob_xpos_x_test_temp = self.Glob_xpos_x_top[where_test_x,j]
				Glob_zpos_x_test_temp = self.Glob_zpos_x_top[where_test_x,j]
				Glob_energy_dep_x_test_temp = self.Glob_energy_dep_x_top[where_test_x,j]
				Glob_pair_flag_x_test_temp = self.Glob_pair_flag_x_top[where_test_x,j]

			where_test_y = np.where(self.Glob_energy_dep_y_top[:,j] > 0.)
			where_test_y = where_test_y[0]
				
			if len(where_test_y) != 0:
				Glob_vol_id_y_test_temp = self.Glob_vol_id_y_top[where_test_y,j]
				Glob_moth_id_y_test_temp = self.Glob_moth_id_y_top[where_test_y,j]
				Glob_Strip_id_y_test_temp = self.Glob_Strip_id_y_top[where_test_y,j]
				Glob_Si_id_y_test_temp = self.Glob_Si_id_y_top[where_test_y,j]
				Glob_tray_id_y_test_temp = self.Glob_tray_id_y_top[where_test_y,j]
				Glob_plane_id_y_test_temp = self.Glob_plane_id_y_top[where_test_y,j]
				Glob_ypos_y_test_temp = self.Glob_ypos_y_top[where_test_y,j]
				Glob_zpos_y_test_temp = self.Glob_zpos_y_top[where_test_y,j]
				Glob_energy_dep_y_test_temp = self.Glob_energy_dep_y_top[where_test_y,j]
				Glob_pair_flag_y_test_temp = self.Glob_pair_flag_y_top[where_test_y,j]
					
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
			while intray < len(Glob_tray_id_test_temp):

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

				intray_max = max(where_tray_eq)
				intray = intray_max + 1
	 	 
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
				event_id_temp_old = self.event_array[j]
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

		self.Glob_event_id_test = Glob_event_id_test 
		self.Glob_vol_id_test = Glob_vol_id_test 
		self.Glob_moth_id_test = Glob_moth_id_test 
		self.Glob_Strip_id_test = Glob_Strip_id_test 
		self.Glob_Si_id_test = Glob_Si_id_test 
		self.Glob_tray_id_test = Glob_tray_id_test 
		self.Glob_plane_id_test = Glob_plane_id_test 
		self.Glob_pos_test = Glob_pos_test 
		self.Glob_zpos_test = Glob_zpos_test 
		self.Glob_energy_dep_test = Glob_energy_dep_test 
		self.Glob_pair_flag_test = Glob_pair_flag_test


########## LO.5 CLUSTER #############

	# IT MAKES THE CLUSTER BY THE BARICENTER METHOD #

    def L05_cluster_x(self):
			
		Glob_event_id_x_top_cluster = []
		Glob_Si_id_x_top_cluster = []
		Glob_tray_id_x_top_cluster = []
		Glob_plane_id_x_top_cluster = []
		Glob_zpos_x_top_cluster = []
		Glob_energy_dep_x_top_cluster = []
		Glob_xpos_x_top_cluster = []
		Glob_Strip_number_x_top_cluster = []
		Glob_pair_flag_x_top_cluster = []
			

		for k in range(self.N_trig):

			N_start = 0
			j=0
				
			while j < np.size(self.Glob_tray_id_x_top[:,k]):

				############ Funzione di ordinamento 3
				
				sort_ascending_plane_x = np.argsort(self.Glob_plane_id_x_top[:,k], kind='mergesort')				
				
				############ fine funzione di ordinamento	


				Glob_vol_id_x_top_tray = self.Glob_vol_id_x_top[sort_ascending_plane_x, k]
				Glob_moth_id_x_top_tray = self.Glob_moth_id_x_top[sort_ascending_plane_x, k]
				Glob_Strip_id_x_top_tray = self.Glob_Strip_id_x_top[sort_ascending_plane_x, k]
				Glob_Si_id_x_top_tray = self.Glob_Si_id_x_top[sort_ascending_plane_x, k]
				Glob_tray_id_x_top_tray = self.Glob_tray_id_x_top[sort_ascending_plane_x, k]
				Glob_plane_id_x_top_tray = self.Glob_plane_id_x_top[sort_ascending_plane_x, k]
				Glob_xpos_x_top_tray = self.Glob_xpos_x_top[sort_ascending_plane_x, k]
				Glob_zpos_x_top_tray = self.Glob_zpos_x_top[sort_ascending_plane_x, k]
				Glob_energy_dep_x_top_tray = self.Glob_energy_dep_x_top[sort_ascending_plane_x, k]
				Glob_pair_flag_x_top_tray = self.Glob_pair_flag_x_top[sort_ascending_plane_x, k]

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


				j_max = max(where_tray_eq_x_top)
				j = j_max + 1


		Glob_event_id_x_top_cluster = np.array((Glob_event_id_x_top_cluster), dtype=np.int64)
		Glob_Si_id_x_top_cluster = np.array((Glob_Si_id_x_top_cluster), dtype=np.int64)			
		Glob_tray_id_x_top_cluster = np.array((Glob_tray_id_x_top_cluster), dtype=np.int64)
		Glob_plane_id_x_top_cluster = np.array((Glob_plane_id_x_top_cluster), dtype=np.int64)
		Glob_zpos_x_top_cluster = np.array(Glob_zpos_x_top_cluster)
		Glob_energy_dep_x_top_cluster = np.array(Glob_energy_dep_x_top_cluster)
		Glob_xpos_x_top_cluster = np.array(Glob_xpos_x_top_cluster)
		Glob_Strip_number_x_top_cluster = np.array((Glob_Strip_number_x_top_cluster), dtype=np.int64)
		Glob_pair_flag_x_top_cluster = np.array(Glob_pair_flag_x_top_cluster)

		self.Glob_event_id_x_top_cluster = Glob_event_id_x_top_cluster 
		self.Glob_Si_id_x_top_cluster = Glob_Si_id_x_top_cluster 
		self.Glob_tray_id_x_top_cluster = Glob_tray_id_x_top_cluster 
		self.Glob_plane_id_x_top_cluster = Glob_plane_id_x_top_cluster 
		self.Glob_zpos_x_top_cluster = Glob_zpos_x_top_cluster 
		self.Glob_energy_dep_x_top_cluster = Glob_energy_dep_x_top_cluster 
		self.Glob_xpos_x_top_cluster = Glob_xpos_x_top_cluster 
		self.Glob_Strip_number_x_top_cluster = Glob_Strip_number_x_top_cluster 
		self.Glob_pair_flag_x_top_cluster = Glob_pair_flag_x_top_cluster


    def L05_cluster_y(self):
	


		Glob_event_id_y_top_cluster = []
		Glob_Si_id_y_top_cluster = []
		Glob_tray_id_y_top_cluster = []
		Glob_plane_id_y_top_cluster = []
		Glob_zpos_y_top_cluster = []
		Glob_energy_dep_y_top_cluster = []
		Glob_ypos_y_top_cluster = []
		Glob_Strip_number_y_top_cluster = []
		Glob_pair_flag_y_top_cluster = []

			
		for k in range(self.N_trig):

			N_start = 0
			j=0

			while j < np.size(self.Glob_tray_id_y_top[:,k]):


				############ Funzione di ordinamento 5
				
				sort_ascending_plane_y = np.argsort(self.Glob_plane_id_y_top[:,k], kind='mergesort')									
				
				############ fine funzione di ordinamento	


				Glob_vol_id_y_top_tray = self.Glob_vol_id_y_top[sort_ascending_plane_y, k]
				Glob_moth_id_y_top_tray = self.Glob_moth_id_y_top[sort_ascending_plane_y, k]
				Glob_Strip_id_y_top_tray = self.Glob_Strip_id_y_top[sort_ascending_plane_y, k]
				Glob_Si_id_y_top_tray = self.Glob_Si_id_y_top[sort_ascending_plane_y, k]
				Glob_tray_id_y_top_tray = self.Glob_tray_id_y_top[sort_ascending_plane_y, k]
				Glob_plane_id_y_top_tray = self.Glob_plane_id_y_top[sort_ascending_plane_y, k]
				Glob_ypos_y_top_tray = self.Glob_ypos_y_top[sort_ascending_plane_y, k]
				Glob_zpos_y_top_tray = self.Glob_zpos_y_top[sort_ascending_plane_y, k]
				Glob_energy_dep_y_top_tray = self.Glob_energy_dep_y_top[sort_ascending_plane_y, k]
				Glob_pair_flag_y_top_tray = self.Glob_pair_flag_y_top[sort_ascending_plane_y, k]

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
										
										
				j_max = max(where_tray_eq_y_top)
				j = j_max + 1

		Glob_event_id_y_top_cluster = np.array((Glob_event_id_y_top_cluster), dtype=np.int64)
		Glob_Si_id_y_top_cluster = np.array((Glob_Si_id_y_top_cluster), dtype=np.int64)			
		Glob_tray_id_y_top_cluster = np.array((Glob_tray_id_y_top_cluster), dtype=np.int64)
		Glob_plane_id_y_top_cluster = np.array((Glob_plane_id_y_top_cluster), dtype=np.int64)
		Glob_zpos_y_top_cluster = np.array(Glob_zpos_y_top_cluster)
		Glob_energy_dep_y_top_cluster = np.array(Glob_energy_dep_y_top_cluster)
		Glob_ypos_y_top_cluster = np.array(Glob_ypos_y_top_cluster)
		Glob_Strip_number_y_top_cluster = np.array((Glob_Strip_number_y_top_cluster), dtype=np.int64)
		Glob_pair_flag_y_top_cluster = np.array(Glob_pair_flag_y_top_cluster)

		self.Glob_event_id_y_top_cluster = Glob_event_id_y_top_cluster 
		self.Glob_Si_id_y_top_cluster = Glob_Si_id_y_top_cluster 
		self.Glob_tray_id_y_top_cluster = Glob_tray_id_y_top_cluster 
		self.Glob_plane_id_y_top_cluster = Glob_plane_id_y_top_cluster
		self.Glob_zpos_y_top_cluster = Glob_zpos_y_top_cluster 
		self.Glob_energy_dep_y_top_cluster = Glob_energy_dep_y_top_cluster 
		self.Glob_ypos_y_top_cluster = Glob_ypos_y_top_cluster 
		self.Glob_Strip_number_y_top_cluster = Glob_Strip_number_y_top_cluster 
		self.Glob_pair_flag_y_top_cluster = Glob_pair_flag_y_top_cluster
		
		
###### L05 MERGING ##############

	# IT PUTS THE X AND Y VIEW OF THE CLUSTER IN THE SAME ARRAY FOR THE L0.5 OUTPUT #	
	
    def merging(self):

		Glob_event_id_cluster = []
		Glob_Si_id_cluster = []
		Glob_tray_id_cluster = []
		Glob_plane_id_cluster = []
		Glob_pos_cluster = []
		Glob_zpos_cluster = []
		Glob_energy_dep_cluster = []
		Glob_Strip_number_cluster = []
		Glob_pair_flag_cluster = []
			
			
		for j in range(self.N_trig):
			Glob_Strip_number_cluster_temp = []
			Glob_Si_id_cluster_temp = []
			Glob_tray_id_cluster_temp = []
			Glob_plane_id_cluster_temp = []
			Glob_pos_cluster_temp = []
			Glob_zpos_cluster_temp = []
			Glob_energy_dep_cluster_temp = []
			Glob_pair_flag_cluster_temp = []

			where_cluster_x_top = np.where(self.Glob_event_id_x_top_cluster == j)
			where_cluster_x_top = where_cluster_x_top[0]

			if len(where_cluster_x_top) != 0:
					
				Glob_Strip_number_cluster_temp_old = self.Glob_Strip_number_x_top_cluster[where_cluster_x_top]
				Glob_Si_id_cluster_temp_old = self.Glob_Si_id_x_top_cluster[where_cluster_x_top]
				Glob_tray_id_cluster_temp_old = self.Glob_tray_id_x_top_cluster[where_cluster_x_top]
				Glob_plane_id_cluster_temp_old = self.Glob_plane_id_x_top_cluster[where_cluster_x_top]
				Glob_pos_cluster_temp_old = self.Glob_xpos_x_top_cluster[where_cluster_x_top]
				Glob_zpos_cluster_temp_old = self.Glob_zpos_x_top_cluster[where_cluster_x_top]
				Glob_energy_dep_cluster_temp_old = self.Glob_energy_dep_x_top_cluster[where_cluster_x_top]
				Glob_pair_flag_cluster_temp_old = self.Glob_pair_flag_x_top_cluster[where_cluster_x_top]

				Glob_Strip_number_cluster_temp = np.append(Glob_Strip_number_cluster_temp, Glob_Strip_number_cluster_temp_old)
				Glob_Si_id_cluster_temp = np.append(Glob_Si_id_cluster_temp, Glob_Si_id_cluster_temp_old)
				Glob_tray_id_cluster_temp = np.append(Glob_tray_id_cluster_temp, Glob_tray_id_cluster_temp_old)
				Glob_plane_id_cluster_temp = np.append(Glob_plane_id_cluster_temp, Glob_plane_id_cluster_temp_old)
				Glob_pos_cluster_temp = np.append(Glob_pos_cluster_temp, Glob_pos_cluster_temp_old)
				Glob_zpos_cluster_temp = np.append(Glob_zpos_cluster_temp, Glob_zpos_cluster_temp_old)
				Glob_energy_dep_cluster_temp = np.append(Glob_energy_dep_cluster_temp, Glob_energy_dep_cluster_temp_old)
				Glob_pair_flag_cluster_temp = np.append(Glob_pair_flag_cluster_temp, Glob_pair_flag_cluster_temp_old)


				where_cluster_y_top = np.where(self.Glob_event_id_y_top_cluster == j)
				where_cluster_y_top = where_cluster_y_top[0]

				if len(where_cluster_y_top) != 0:

					Glob_Strip_number_cluster_temp_old = self.Glob_Strip_number_y_top_cluster[where_cluster_y_top]
					Glob_Si_id_cluster_temp_old = self.Glob_Si_id_y_top_cluster[where_cluster_y_top]
					Glob_tray_id_cluster_temp_old = self.Glob_tray_id_y_top_cluster[where_cluster_y_top]
					Glob_plane_id_cluster_temp_old = self.Glob_plane_id_y_top_cluster[where_cluster_y_top]
					Glob_pos_cluster_temp_old = self.Glob_ypos_y_top_cluster[where_cluster_y_top]
					Glob_zpos_cluster_temp_old = self.Glob_zpos_y_top_cluster[where_cluster_y_top]
					Glob_energy_dep_cluster_temp_old = self.Glob_energy_dep_y_top_cluster[where_cluster_y_top]
					Glob_pair_flag_cluster_temp_old = self.Glob_pair_flag_y_top_cluster[where_cluster_y_top]

					Glob_Strip_number_cluster_temp = np.append(Glob_Strip_number_cluster_temp, Glob_Strip_number_cluster_temp_old)
					Glob_Si_id_cluster_temp = np.append(Glob_Si_id_cluster_temp, Glob_Si_id_cluster_temp_old)
					Glob_tray_id_cluster_temp = np.append(Glob_tray_id_cluster_temp, Glob_tray_id_cluster_temp_old)
					Glob_plane_id_cluster_temp = np.append(Glob_plane_id_cluster_temp, Glob_plane_id_cluster_temp_old)
					Glob_pos_cluster_temp = np.append(Glob_pos_cluster_temp, Glob_pos_cluster_temp_old)
					Glob_zpos_cluster_temp = np.append(Glob_zpos_cluster_temp, Glob_zpos_cluster_temp_old)
					Glob_energy_dep_cluster_temp = np.append(Glob_energy_dep_cluster_temp, Glob_energy_dep_cluster_temp_old)
					Glob_pair_flag_cluster_temp = np.append(Glob_pair_flag_cluster_temp, Glob_pair_flag_cluster_temp_old)

			else:
					
				where_cluster_y_top = np.where(self.Glob_event_id_y_top_cluster == j)
				where_cluster_y_top = where_cluster_y_top[0]

				if len(where_cluster_y_top) != 0:

					Glob_Strip_number_cluster_temp_old = self.Glob_Strip_number_y_top_cluster[where_cluster_y_top]
					Glob_Si_id_cluster_temp_old = self.Glob_Si_id_y_top_cluster[where_cluster_y_top]
					Glob_tray_id_cluster_temp_old = self.Glob_tray_id_y_top_cluster[where_cluster_y_top]
					Glob_plane_id_cluster_temp_old = self.Glob_plane_id_y_top_cluster[where_cluster_y_top]
					Glob_pos_cluster_temp_old = self.Glob_ypos_y_top_cluster[where_cluster_y_top]
					Glob_zpos_cluster_temp_old = self.Glob_zpos_y_top_cluster[where_cluster_y_top]
					Glob_energy_dep_cluster_temp_old = self.Glob_energy_dep_y_top_cluster[where_cluster_y_top]
					Glob_pair_flag_cluster_temp_old = self.Glob_pair_flag_y_top_cluster[where_cluster_y_top]

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
			while intray < len(Glob_tray_id_cluster_temp):
				
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
					

				intray_max = max(where_tray_eq)
				intray = intray_max + 1

			event_id_temp = np.zeros(len(Si_id_temp), dtype=np.int64)
			for k in range(len(Si_id_temp)):
				event_id_temp[k] = self.event_array[j]

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

		self.Glob_event_id_cluster = Glob_event_id_cluster 
		self.Glob_Si_id_cluster = Glob_Si_id_cluster 
		self.Glob_tray_id_cluster = Glob_tray_id_cluster 
		self.Glob_plane_id_cluster = Glob_plane_id_cluster 
		self.Glob_pos_cluster = Glob_pos_cluster 
		self.Glob_zpos_cluster = Glob_zpos_cluster 
		self.Glob_energy_dep_cluster = Glob_energy_dep_cluster
		self.Glob_Strip_number_cluster = Glob_Strip_number_cluster 
		self.Glob_pair_flag_cluster = Glob_pair_flag_cluster


########## ENERGY RESOLUTION ##########
		#It computes the discards between the position of the cluster and the hit in every plane in order to show the energy resolution of the tracker
		
	##Cluster.dat
    def resolution(self, N_in, part_type, ang_type, theta_type, phi_type, ifile):
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('                 Tracker Resolution   ')
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		
		if os.path.exists(self.outdir+'/'+self.sim_tag+'_RESOLUTION_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_RESOLUTION_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data = open(self.outdir+'/'+self.sim_tag+'_RESOLUTION_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data = open(self.outdir+'/'+self.sim_tag+'_RESOLUTION_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		#self.outdir+'/'+self.sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'
		#file cluster
		event_x_list = []
		event_y_list = []
		plane_x_list = []
		plane_y_list = []
		pos_x_list = []
		pos_y_list = []
		si_x_list = []
		si_y_list = []
		
		j = 0
		while j < len(self.Glob_event_id_cluster):
			where_event_eq = np.where(self.Glob_event_id_cluster == self.Glob_event_id_cluster[j]) 
			where_event_eq = where_event_eq[0]

			event = self.Glob_event_id_cluster[where_event_eq]
			plane = self.Glob_plane_id_cluster[where_event_eq] 
			#zpos = self.Glob_zpos_cluster[where_event_eq]
			si = self.Glob_Si_id_cluster[where_event_eq]
			pos = self.Glob_pos_cluster[where_event_eq]
			#en_dep = self.Glob_energy_dep_cluster[where_event_eq]
			n_strips = self.Glob_Strip_number_cluster[where_event_eq]
			#pair_flag = self.Glob_pair_flag_cluster[where_event_eq]

			#------------------
			# X VIEW
			#------------------

			where_x = np.where((si == 0) & (n_strips > 1)) #elementi formanti il cluster
			where_x = where_x[0]

			event_x_list.append(event[where_x[0]])
			plane_x_list.append(plane[where_x[0]])
			pos_x_list.append(pos[where_x[0]])
			si_x_list.append(si[where_x[0]])
					
					
			data.write('{:d}\t'.format(event[where_x[0]]))
			data.write('{:d}\t'.format(si[where_x[0]]))
			data.write('{:d}\t'.format(plane[where_x[0]]))
			data.write('{:f}\t'.format(pos[where_x[0]]))
			data.write('{:d}\n'.format(n_strips[where_x[0]]))
			
			#------------------
			# Y VIEW
			#------------------

#			where_y = np.where((si == 1) & (n_strips > 1))
#			where_y = where_y[0]

					
#			event_y_list.append(event[where_y[0]])
#			plane_y_list.append(plane[where_x[0]])
#			pos_y_list.append(pos[where_x[0]])
#			si_y_list.append(si[where_x[0]])
				
#			data.write('{:d}\t'.format(event[where_y[0]]))
#			data.write('{:d}\t'.format(si[where_y[0]]))
#			data.write('{:d}\t'.format(plane[where_y[0]]))
#			data.write('{:f}\t'.format(pos[where_y[0]]))
#			data.write('{:d}\n'.format(n_strips[where_y[0]]))
							
											
			j_max = max(where_event_eq)
			j = j_max + 1
			
		data.close()
		
			
		#file AA strip.dat
		
		if os.path.exists(self.outdir+'/'+self.sim_tag+'_RESOLUTION_STRIP_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_RESOLUTION_STRIP_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data1 = open(self.outdir+'/'+self.sim_tag+'_RESOLUTION_STRIP_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data1 = open(self.outdir+'/'+self.sim_tag+'_RESOLUTION_STRIP_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		event_0_x_list = []
		event_0_y_list = []
		plane_0_x_list = []
		plane_0_y_list = []
		pos_0_x_list = []
		pos_0_y_list = []
		si_0_x_list = []
		si_0_y_list = []
		zpos_0_x_list = []
		zpos_0_y_list = []
		
		j = 0
		while j < len(self.Glob_event_id_test):
			where_event_eq_0 = np.where(self.Glob_event_id_test == self.Glob_event_id_test[j]) 
			where_event_eq_0 = where_event_eq_0[0]

			event_0 = self.Glob_event_id_test[where_event_eq_0]
			plane_0 = self.Glob_plane_id_test[where_event_eq_0] 
			zpos_0 = self.Glob_zpos_test[where_event_eq_0]
			si_0 = self.Glob_Si_id_test[where_event_eq_0]
			pos_0 = self.Glob_pos_test[where_event_eq_0]
			#en_dep = self.Glob_energy_dep_cluster[where_event_eq]
			strip_id_0 = self.Glob_Strip_id_test[where_event_eq_0]
			#pair_flag = self.Glob_pair_flag_cluster[where_event_eq]

			#------------------
			# X VIEW
			#------------------
			where_x_0 = np.where(si_0 == 0)
			where_x_0 = where_x_0[0]
			
			r = 0
			if len(where_x_0) != 0:
				while r < len(where_x_0):
					if ((strip_id_0[where_x_0[r]] == strip_id_0[where_x_0[r]-1]-1) & (si_0[where_x_0[r]-1] == 0)):
						event_0_x = event_0[where_x_0[r]-1]
						plane_0_x = plane_0[where_x_0[r]-1]
						pos_0_x = pos_0[where_x_0[r]-1]
						si_0_x = si_0[where_x_0[r]-1]
						zpos_0_x = zpos_0[where_x_0[r]-1]

						event_0_x_list.append(event_0_x)
						plane_0_x_list.append(plane_0_x)
						pos_0_x_list.append(pos_0_x)
						si_0_x_list.append(si_0_x)
						zpos_0_x_list.append(zpos_0_x)
						
						data1.write('{:d}\t'.format(event_0[where_x_0[r]-1]))
						data1.write('{:d}\t'.format(si_0[where_x_0[r]-1]))
						data1.write('{:d}\t'.format(plane_0[where_x_0[r]-1]))
						data1.write('{:f}\t'.format(pos_0[where_x_0[r]-1]))
						data1.write('{:d}\n'.format(strip_id_0[where_x_0[r]-1])) 
							
						break
					
					r = r + 1

			#------------------
			# Y VIEW
			#------------------
#			where_y_0 = np.where(si_0 == 1)
#			where_y_0 = where_y_0[0]
			
#			r = 0
#			if len(where_y_0) != 0:
#				while r < len(where_y_0):
#					if ((strip_id_0[where_y_0[r]] == strip_id_0[where_y_0[r]-1]-1) & (si_0[where_y_0[r]-1] == 1)):
#						event_0_y = event_0[where_y_0[r]-1]
#						plane_0_y = plane_0[where_y_0[r]-1]
#						pos_0_y = pos_0[where_y_0[r]-1]
#						si_0_y = si_0[where_y_0[r]-1]
#						zpos_0_y = zpos_0[where_y_0[r]-1]
				
#						event_0_y_list.append(event_0_y)
#						plane_0_y_list.append(plane_0_y)
#						pos_0_y_list.append(pos_0_y)
#						si_0_y_list.append(si_0_y)
#						zpos_0_y_list.append(zpos_0_y)
						
#						data1.write('{:d}\t'.format(event_0[where_y_0[r]-1]))
#						data1.write('{:d}\t'.format(si_0[where_y_0[r]-1]))
#						data1.write('{:d}\t'.format(plane_0[where_y_0[r]-1]))
#						data1.write('{:f}\t'.format(pos_0[where_y_0[r]-1]))
#						data1.write('{:d}\n'.format(strip_id_0[where_y_0[r]-1])) 
					
#						break
					
#					r = r + 1
				
			j_max = max(where_event_eq_0)
			j = j_max + 1
		
		data1.close()
		
		print('cluster',len(pos_x_list))
#		print('cluster_y',len(pos_y_list))
		print('strip',len(pos_0_x_list))
#		print('strip_y',len(pos_0_y_list))
		
		event_x_list = np.array(event_x_list)
#		event_y_list = np.array(event_y_list)
		plane_x_list = np.array(plane_x_list)
#		plane_y_list = np.array(plane_y_list)
		pos_x_list = np.array(pos_x_list)
#		pos_y_list = np.array(pos_y_list)
		si_x_list = np.array(si_x_list)
#		si_y_list = np.array(si_y_list)
		
		event_0_x_list = np.array(event_0_x_list)
#		event_0_y_list = np.array(event_0_y_list)
		plane_0_x_list = np.array(plane_0_x_list)
#		plane_0_y_list = np.array(plane_0_y_list)
		pos_0_x_list = np.array(pos_0_x_list)
#		pos_0_y_list = np.array(pos_0_y_list)
		si_0_x_list = np.array(si_0_x_list)
#		si_0_y_list = np.array(si_0_y_list)

		#### Difference X #######

		if os.path.exists(self.outdir+'/'+self.sim_tag+'_DIFFERENCE_RESOLUTION_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat'):
			os.remove(self.outdir+'/'+self.sim_tag+'_DIFFERENCE_RESOLUTION_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat')
			data_diff_x = open(self.outdir+'/'+self.sim_tag+'_DIFFERENCE_RESOLUTION_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		else:
			data_diff_x = open(self.outdir+'/'+self.sim_tag+'_DIFFERENCE_RESOLUTION_'+str(N_in)+part_type+'_'+self.sname+'_'+self.ene_dis+'_'+ang_type+'_'+self.ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+self.pol_string+str(ifile)+'.dat', 'w')
		
		diff_x_list = []
		t = 0 
		while t < len(pos_x_list):
			diff_x = pos_0_x_list[t]-pos_x_list[t]
			diff_x_list.append(diff_x)
			
			data_diff_x.write('{:d}\t'.format(event_x_list[t]))
			data_diff_x.write('{:d}\t'.format(plane_x_list[t]))
			data_diff_x.write('{:f}\n'.format(diff_x))

			t = t + 1
		diff_x_list_abs = np.absolute(diff_x_list)
		mvx = (sum(diff_x_list_abs)/len(pos_x_list)*1000.)  #mm ---> micron
		print('mean_value',mvx,'micron')
		
		data_diff_x.close()
		

		#### Difference Y #######

#		if os.path.exists('difference_y.dat'):
#			os.remove('difference_y.dat')
#			data_diff_y = open('difference_y.dat', 'w')
#		else:
#			data_diff_y = open('difference_y.dat', 'w')

#		diff_y_list = []
#		t = 0 
#		while t < len(pos_y_list):
#			diff_y = pos_0_y_list[t]-pos_y_list[t]
#			diff_y_list.append(diff_y)

#			data_diff_y.write('{:d}\t'.format(event_y_list[t]))
#			data_diff_y.write('{:d}\t'.format(plane_y_list[t]))
#			data_diff_y.write('{:f}\n'.format(diff_y))

#			t = t + 1
#		diff_y_list_abs = np.absolute(diff_y_list)
#		mvy = (sum(diff_y_list_abs)/len(pos_y_list)*1000.)  #mm ---> micron
#		print('mean_value_y',mvy,'micron')
		
#		data_diff_y.close()	

				
########## CALORIMETER ###########

	# IT FLAGS THE PRIMARY EVENT, SEPARATING COMPTON FROM PAIR #
	# IT SUMS THE ENERGY OF THE SAME EVENT FOR EVERY VOLUME #
	# IT FLAGS THE SUMMED EVENT IF ONE HITS ARE IN THE SAME VOLUME OR IT FLAGS THE SUMMED EVENT IF THE HITS AREN'T IN THE SAME VOLUME  #

    def G4_cal(self):

		bar_ene = self.energy_dep_cal

		N_trig_cal = 0

		event_id_tot_cal = []
		vol_id_tot_cal = []
		bar_id_tot = []
		moth_id_tot_cal = []
		bar_ene_tot = []
		pair_flag_tot_cal = []
		
		j=0
		while j < len(self.event_id_cal):
			
			where_event_eq = np.where(self.event_id_cal == self.event_id_cal[j])
			where_event_eq = where_event_eq[0]

			N_trig_cal = N_trig_cal + 1
			
			vol_id_temp_cal = self.vol_id_cal[where_event_eq]
			moth_id_temp_cal = self.moth_id_cal[where_event_eq]
			bar_ene_temp = bar_ene[where_event_eq]
			trk_id_temp_cal = self.trk_id_cal[where_event_eq]
			child_id_temp_cal = self.child_id_cal[where_event_eq]
			proc_id_temp_cal = self.proc_id_cal[where_event_eq]
			gtime_temp_cal = self.gtime_ent_cal[where_event_eq]

			where_trk_event_cal = np.where(self.event_id_tot_tr_raw == self.event_id_cal[j])
			where_trk_event_cal = where_trk_event_cal[0]
			if (where_trk_event_cal.size):
				cal_event_flag = self.pair_flag_tot_tr_raw[where_trk_event_cal]
			else:
				cal_event_flag = 0
			
			"""
			where_pair = np.where((child_id_temp_cal == 1) & (proc_id_temp_cal == 7) & (trk_id_temp_cal <= 3))
			where_pair = where_pair[0]
			where_compton = np.where((child_id_temp_cal == 1) & (proc_id_temp_cal == 3) & (trk_id_temp_cal <= 2))
			where_compton = where_compton[0]
			where_ray = np.where((child_id_temp_cal == 1) & (proc_id_temp_cal == 9) & (trk_id_temp_cal == 1))
			where_ray = where_ray[0]
			
			ispair = 0
			iscompton = 0
			isray = 0
			isother = 0
			gtime_pair = [10**9]
			gtime_compton = [10**9]
			gtime_ray = [10**9]
			
			if len(where_pair) != 0:
				gtime_pair = gtime_temp[where_pair]
				#ispair = 1
			else:
				gtime_pair = [10**9]
				#ispair = 0

			if len(where_compton) != 0:
				gtime_compton = gtime_temp[where_compton]
				#iscompton = 1
			else:
				gtime_compton = [10**9]
				#iscompton = 0

			if len(where_ray) != 0:
				gtime_ray = gtime_temp[where_ray]
				#isray = 1
			else:
				gtime_ray = [10**9]
				#isray = 0
			
			if ((gtime_pair[0] == 10**9) and (gtime_compton[0] == 10**9) and (gtime_ray[0] == 10**9)):
				isother = 1
			else:
				gtime_array = [gtime_pair[0], gtime_compton[0], gtime_ray[0]]
				proc_index = np.argmin(gtime_array)
				if proc_index == 0: ispair = 1
				if proc_index == 1: iscompton = 1
				if proc_index == 2: isray = 1
			"""

			r = 0
			while 1:
				ispair_vol = 0
				iscompton_vol = 0
				isray_vol = 0
				isother_vol = 0
				isprimary_vol = 0
	
				where_vol_eq = np.where((vol_id_temp_cal == vol_id_temp_cal[r]) & (moth_id_temp_cal == moth_id_temp_cal[r]))										
				where_vol_eq = where_vol_eq[0]	

				where_other_vol = np.where((vol_id_temp_cal != vol_id_temp_cal[r]) | (moth_id_temp_cal != moth_id_temp_cal[r]))
				where_other_vol = where_other_vol[0]
				
				bar_ene_tot_temp = np.sum(bar_ene_temp[where_vol_eq])
				
				if bar_ene_tot_temp >= self.E_th_cal:
					event_id_tot_cal_old = self.event_id_cal[j]
					vol_id_tot_cal_old = vol_id_temp_cal[r]
					bar_id_tot_old = vol_id_temp_cal[r] - self.cal_vol_start
					moth_id_tot_cal_old = moth_id_temp_cal[r]
					bar_ene_tot_old = np.sum(bar_ene_temp[where_vol_eq])

					event_id_tot_cal.append(event_id_tot_cal_old)
					vol_id_tot_cal.append(vol_id_tot_cal_old)
					bar_id_tot.append(bar_id_tot_old )
					moth_id_tot_cal.append(moth_id_tot_cal_old)
					bar_ene_tot.append(bar_ene_tot_old) 					
					
					#Searching for Pair/Compton events of secondary particles generated by the primary
					#if another process is involved, the event is flagged as 0
					#if one of hits in the same volume is a pair the summed event is flagged as 1
					#if one of hits in the same volume is a compton the summed event is flagged as 2
					#if one of hits in the same volume is a rayleigh the summed event is flagged as 3

					all_child = child_id_temp_cal[where_vol_eq]
					all_trk = trk_id_temp_cal[where_vol_eq]
					all_proc = proc_id_temp_cal[where_vol_eq]
					all_gtime = gtime_temp_cal[where_vol_eq]
					
					ispair_vol = 0
					iscompton_vol = 0					
					
					where_pair_vol = np.where((all_child == 1) & (all_proc == 7) & (all_trk <= 3))
					where_pair_vol = where_pair_vol[0]
					
					where_compton_vol = np.where(((all_child == 1) & (all_proc == 3) & (all_trk == 2)) | ((all_child == 0) & (all_trk == 1)))
					where_compton_vol = where_compton_vol[0]
					
					#where_ray_vol = np.where((all_child == 1) & (all_proc == 9) & (all_trk == 1))
					#where_ray_vol = where_ray_vol[0]
					where_pair_in_trk = np.where(cal_event_flag == 1)
					where_pair_in_trk = where_pair_in_trk[0]
					if len(where_pair_vol) != 0 & (len(where_pair_in_trk) != 0):
						pair_flag_tot_old = 1
						ispair_vol = 1
						
						pair_flag_tot_cal.append(pair_flag_tot_old)	

					where_compton_in_trk = np.where(cal_event_flag == 2)
					where_compton_in_trk = where_compton_in_trk[0]						
					if len(where_compton_vol) != 0 & (len(where_compton_in_trk) != 0):
						if ispair_vol == 0:
							pair_flag_tot_old = 2
							iscompton_vol = 1
						
							pair_flag_tot_cal.append(pair_flag_tot_old)	
						
					if ((ispair_vol == 0) & (iscompton_vol == 0)):
						if ((len(where_pair_vol) == 0) & (len(where_compton_vol) == 0)):
							where_ray_in_trk = np.where(cal_event_flag == 3)
							where_ray_in_trk = where_ray_in_trk[0]						
							if (len(where_ray_in_trk) != 0):
								pair_flag_tot_old = 3
								isray_vol = 1
							
								pair_flag_tot_cal.append(pair_flag_tot_old)	

							else:
								pair_flag_tot_old = 0
								pair_flag_tot_cal.append(pair_flag_tot_old)

						else:
							if ((len(where_pair_vol) != 0 & (len(where_pair_in_trk) == 0)) | (len(where_compton_vol) != 0 & (len(where_compton_in_trk) == 0))):
								pair_flag_tot_old = 0
								pair_flag_tot_cal.append(pair_flag_tot_old)


				if len(where_other_vol) != 0:
					vol_id_temp_cal = vol_id_temp_cal[where_other_vol]
					moth_id_temp_cal = moth_id_temp_cal[where_other_vol]
					bar_ene_temp = bar_ene_temp[where_other_vol]
					trk_id_temp_cal = trk_id_temp_cal[where_other_vol]
					child_id_temp_cal = child_id_temp_cal[where_other_vol]
					proc_id_temp_cal = proc_id_temp_cal[where_other_vol]
					gtime_temp_cal = gtime_temp_cal[where_other_vol]
				
				else:
					break
						

			j_max = max(where_event_eq)
			j = j_max + 1
							
		event_id_tot_cal = np.array((event_id_tot_cal), dtype = np.int64)
		bar_id_tot = np.array((bar_id_tot), dtype = np.int64)
		bar_ene_tot = np.array(bar_ene_tot)
		pair_flag_tot_cal = np.array(pair_flag_tot_cal)
		
		self.event_id_tot_cal = event_id_tot_cal 
		self.bar_id_tot = bar_id_tot 
		self.bar_ene_tot = bar_ene_tot
		self.pair_flag_tot_cal = pair_flag_tot_cal


	######## COMPTON/PAIR CAL EVENTS ######
	
	# IT SEPARATES THE COMPTON EVENTS FROM THE PAIR PRODUCTION ONES
	
    def compton_cal(self):

		where_compton_cal = np.where(self.pair_flag_tot_cal == 2)
		where_compton_cal = where_compton_cal[0]

		if where_compton_cal.size:
			event_id_tot_cal_compton = self.event_id_tot_cal[where_compton_cal]
			bar_id_tot_compton = self.bar_id_tot[where_compton_cal]
			bar_ene_tot_compton = self.bar_ene_tot[where_compton_cal]
			pair_flag_tot_cal_compton = self.pair_flag_tot_cal[where_compton_cal]

			self.event_id_tot_cal_compton = event_id_tot_cal_compton
			self.bar_id_tot_compton = bar_id_tot_compton
			self.bar_ene_tot_compton = bar_ene_tot_compton
			self.pair_flag_tot_cal_compton = pair_flag_tot_cal_compton
			
	
		else:
			print('No Compton CAL events found')
	
		self.where_compton_cal = where_compton_cal 

		
    def pair_cal(self):

		where_pair_cal = np.where(self.pair_flag_tot_cal == 1)
		where_pair_cal = where_pair_cal[0]

		if where_pair_cal.size:
			event_id_tot_cal_pair = self.event_id_tot_cal[where_pair_cal]
			bar_id_tot_pair = self.bar_id_tot[where_pair_cal]
			bar_ene_tot_pair = self.bar_ene_tot[where_pair_cal]
			pair_flag_tot_cal_pair = self.pair_flag_tot_cal[where_pair_cal]

			self.event_id_tot_cal_pair = event_id_tot_cal_pair
			self.bar_id_tot_pair = bar_id_tot_pair
			self.bar_ene_tot_pair = bar_ene_tot_pair
			self.pair_flag_tot_cal_pair = pair_flag_tot_cal_pair
		


		else:
			print('No Pair CAL events found')

		self.where_pair_cal = where_pair_cal

    def ray_cal(self):

		where_ray_cal = np.where(self.pair_flag_tot_cal == 3)
		where_ray_cal = where_ray_cal[0]

		if where_ray_cal.size:
			event_id_tot_cal_ray = self.event_id_tot_cal[where_ray_cal]
			bar_id_tot_ray = self.bar_id_tot[where_ray_cal]
			bar_ene_tot_ray = self.bar_ene_tot[where_ray_cal]
			pair_flag_tot_cal_ray = self.pair_flag_tot_cal[where_ray_cal]

			self.event_id_tot_cal_ray = event_id_tot_cal_ray
			self.bar_id_tot_ray = bar_id_tot_ray
			self.bar_ene_tot_ray = bar_ene_tot_ray
			self.pair_flag_tot_cal_ray = pair_flag_tot_cal_ray
		


		else:
			print('No Rayleigh CAL events found')

		self.where_ray_cal = where_ray_cal
		
	################################

	# IT SUMS THE ENERGY FOR EVERY EVENT #  

    def cal_sum(self):

		event_id_tot_cal_sum = []
		bar_ene_tot_sum = []


		j=0
		while j < len(self.event_id_tot_cal):
			
			where_event_eq = np.where(self.event_id_tot_cal == self.event_id_tot_cal[j])
			where_event_eq = where_event_eq[0]
			
			
			event_id_tot_cal_sum_old = self.event_id_tot_cal[j]
			
			event_id_tot_cal_sum.append(event_id_tot_cal_sum_old)
			
			bar_ene_tot_sum_old = np.sum(self.bar_ene_tot[where_event_eq])
			bar_ene_tot_sum = np.append(bar_ene_tot_sum, bar_ene_tot_sum_old)

			j_max = max(where_event_eq)
			j = j_max + 1

		event_id_tot_cal_sum = np.array((event_id_tot_cal_sum), dtype=np.int64)
		bar_ene_tot_sum = np.array(bar_ene_tot_sum)
	
		self.event_id_tot_cal_sum = event_id_tot_cal_sum 
		self.bar_ene_tot_sum = bar_ene_tot_sum
	
		
###### ANTICOINCIDENCE ######

	# IT FLAGS THE PRIMARY EVENT, SEPARATING COMPTON FROM PAIR #
	# IT SUMS THE ENERGY OF THE SAME EVENT FOR EVERY VOLUME #
	# IT FLAGS THE SUMMED EVENT IF ONE HITS ARE IN THE SAME VOLUME OR IT FLAGS THE SUMMED EVENT IF THE HITS AREN'T IN THE SAME VOLUME  #
	# IT CREATES THE AC PANEL AND SUBPANEL ID

    def AC_analysis(self):

		N_trig_ac = 0


		event_id_tot_ac = []
		vol_id_tot_ac = []
		moth_id_tot_ac = []
		energy_dep_tot_ac = []
		pair_flag_tot_ac = []


		j=0
		while j < len(self.event_id_ac):
		
			where_event_eq = np.where(self.event_id_ac == self.event_id_ac[j])
			where_event_eq = where_event_eq[0]

			N_trig_ac = N_trig_ac + 1

			vol_id_temp_ac = self.vol_id_ac[where_event_eq]
			moth_id_temp_ac = self.moth_id_ac[where_event_eq]
			energy_dep_temp_ac = self.energy_dep_ac[where_event_eq]
			trk_id_temp_ac = self.trk_id_ac[where_event_eq]
			child_id_temp_ac = self.child_id_ac[where_event_eq]
			proc_id_temp_ac = self.proc_id_ac[where_event_eq]
			gtime_temp_ac = self.gtime_ent_ac[where_event_eq]


			where_trk_event_ac = np.where(self.event_id_tot_tr_raw == self.event_id_ac[j])
			where_trk_event_ac = where_trk_event_ac[0]
			if (where_trk_event_ac.size):
				ac_event_flag = self.pair_flag_tot_tr_raw[where_trk_event_ac]
			else:
				ac_event_flag = 0
				
			"""
			where_pair = np.where((child_id_temp_ac == 1) & (proc_id_temp_ac == 7) & (trk_id_temp_ac <= 3))
			where_pair = where_pair[0]
			where_compton = np.where((child_id_temp_ac == 1) & (proc_id_temp_ac == 3) & (trk_id_temp_ac <= 2))
			where_compton = where_compton[0]
			where_ray = np.where((child_id_temp_ac == 1) & (proc_id_temp_ac == 9) & (trk_id_temp_ac == 1))
			where_ray = where_ray[0]
			
			ispair = 0
			iscompton = 0
			isray = 0
			isother = 0
			gtime_pair = [10**9]
			gtime_compton = [10**9]
			gtime_ray = [10**9]
			
			if len(where_pair) != 0:
				gtime_pair = gtime_temp[where_pair]
				#ispair = 1
			else:
				gtime_pair = [10**9]
				#ispair = 0

			if len(where_compton) != 0:
				gtime_compton = gtime_temp[where_compton]
				#iscompton = 1
			else:
				gtime_compton = [10**9]
				#iscompton = 0

			if len(where_ray) != 0:
				gtime_ray = gtime_temp[where_ray]
				#isray = 1
			else:
				gtime_ray = [10**9]
				#isray = 0
			
			if ((gtime_pair[0] == 10**9) and (gtime_compton[0] == 10**9) and (gtime_ray[0] == 10**9)):
				isother = 1
			else:
				gtime_array = [gtime_pair[0], gtime_compton[0], gtime_ray[0]]
				proc_index = np.argmin(gtime_array)
				if proc_index == 0: ispair = 1
				if proc_index == 1: iscompton = 1
				if proc_index == 2: isray = 1
			"""


			r = 0
			while 1:

				ispair_vol = 0
				iscompton_vol = 0
				isray_vol = 0
				isother_vol = 0
				isprimary_vol = 0
	
				where_vol_eq = np.where(vol_id_temp_ac == vol_id_temp_ac[r])
				where_vol_eq = where_vol_eq[0]

				where_other_vol = np.where(vol_id_temp_ac != vol_id_temp_ac[r])
				where_other_vol = where_other_vol[0]
				
				event_id_tot_ac_old = self.event_id_ac[j]
				vol_id_tot_ac_old = vol_id_temp_ac[r]
				moth_id_tot_ac_old = moth_id_temp_ac[r]
				energy_dep_tot_ac_old = np.sum(energy_dep_temp_ac[where_vol_eq])
	
				event_id_tot_ac.append(event_id_tot_ac_old)
				vol_id_tot_ac.append(vol_id_tot_ac_old)
				moth_id_tot_ac.append(moth_id_tot_ac_old)
				energy_dep_tot_ac.append(energy_dep_tot_ac_old)				
				
				#Searching for Pair/Compton events of secondary particles generated by the primary
				#if another process is involved, the event is flagged as 0
				#if one of hits in the same volume is a pair the summed event is flagged as 1
				#if one of hits in the same volume is a compton the summed event is flagged as 2
				#if one of hits in the same volume is a rayleigh the summed event is flagged as 3
				
				all_child = child_id_temp_ac[where_vol_eq]
				all_trk = trk_id_temp_ac[where_vol_eq]
				all_proc = proc_id_temp_ac[where_vol_eq]
				all_gtime = gtime_temp_ac[where_vol_eq]
				
				ispair_vol = 0
				iscompton_vol = 0
				
				where_pair_vol = np.where((all_child == 1) & (all_proc == 7) & (all_trk <= 3))
				where_pair_vol = where_pair_vol[0]
				
				where_compton_vol = np.where(((all_child == 1) & (all_proc == 3) & (all_trk == 2)) | ((all_child == 0) & (all_trk == 1)))
				where_compton_vol = where_compton_vol[0]
				
				#where_ray_vol = np.where((all_child == 1) & (all_proc == 9) & (all_trk == 1))
				#where_ray_vol = where_ray_vol[0]
				
				where_pair_in_trk = np.where(ac_event_flag == 1)
				where_pair_in_trk = where_pair_in_trk[0]
				if len(where_pair_vol) != 0 & (len(where_pair_in_trk) != 0):
					pair_flag_tot_old = 1
					ispair_vol = 1
					
					pair_flag_tot_ac.append(pair_flag_tot_old)	

				where_compton_in_trk = np.where(ac_event_flag == 2)
				where_compton_in_trk = where_compton_in_trk[0]						
				if len(where_compton_vol) != 0 & (len(where_compton_in_trk) != 0):
					if ispair_vol == 0:
						pair_flag_tot_old = 2
						iscompton_vol = 1
					
						pair_flag_tot_ac.append(pair_flag_tot_old)	
		
				if ((ispair_vol == 0) & (iscompton_vol == 0)):
					if ((len(where_pair_vol) == 0) & (len(where_compton_vol) == 0)):
						where_ray_in_trk = np.where(ac_event_flag == 3)
						where_ray_in_trk = where_ray_in_trk[0]						
						if (len(where_ray_in_trk) != 0):
							pair_flag_tot_old = 3
							isray_vol = 1
						
							pair_flag_tot_ac.append(pair_flag_tot_old)	
						else:
							pair_flag_tot_old = 0
							pair_flag_tot_ac.append(pair_flag_tot_old)
					else:
						if ((len(where_pair_vol) != 0 & (len(where_pair_in_trk) == 0)) | (len(where_compton_vol) != 0 & (len(where_compton_in_trk) == 0))):
							pair_flag_tot_old = 0
							pair_flag_tot_ac.append(pair_flag_tot_old)
					

				if len(where_other_vol) != 0:
					vol_id_temp_ac = vol_id_temp_ac[where_other_vol]
					moth_id_temp_ac = moth_id_temp_ac[where_other_vol]
					energy_dep_temp_ac = energy_dep_temp_ac[where_other_vol]
					trk_id_temp_ac = trk_id_temp_ac[where_other_vol]
					child_id_temp_ac = child_id_temp_ac[where_other_vol]
					proc_id_temp_ac = proc_id_temp_ac[where_other_vol]
					gtime_temp_ac = gtime_temp_ac[where_other_vol]
				
				else:
					break

			j_max = max(where_event_eq)
			j = j_max + 1
		
		vol_id_tot_ac = np.array(vol_id_tot_ac)

		# AC panel IDs

		AC_panel = [' ']*(len(vol_id_tot_ac)) 
		AC_subpanel = [0]*(len(vol_id_tot_ac))

		for j in range(len(vol_id_tot_ac)):

			if vol_id_tot_ac[j] >= self.panel_S[0] and vol_id_tot_ac[j] <= self.panel_S[2]:
				AC_panel[j] = 'S'
				
				if vol_id_tot_ac[j] == self.panel_S[0]:
					AC_subpanel[j] = 3
					
				if vol_id_tot_ac[j] == self.panel_S[1]:
					AC_subpanel[j] = 2
				
				if vol_id_tot_ac[j] == self.panel_S[2]:
					AC_subpanel[j] = 1
			
					

			if vol_id_tot_ac[j] >= self.panel_D[0] and vol_id_tot_ac[j] <= self.panel_D[2]:
				AC_panel[j] = 'D'
				
				if vol_id_tot_ac[j] == self.panel_D[0]:
					AC_subpanel[j] = 3
					
				if vol_id_tot_ac[j] == self.panel_D[1]:
					AC_subpanel[j] = 2
				
				if vol_id_tot_ac[j] == self.panel_D[2]:
					AC_subpanel[j] = 1


			if vol_id_tot_ac[j] >= self.panel_F[0] and vol_id_tot_ac[j] <= self.panel_F[2]:
				AC_panel[j] = 'F'
				
				if vol_id_tot_ac[j] == self.panel_F[0]:
					AC_subpanel[j] = 1
					
				if vol_id_tot_ac[j] == self.panel_F[1]:
					AC_subpanel[j] = 2
				
				if vol_id_tot_ac[j] == self.panel_F[2]:
					AC_subpanel[j] = 3


			if vol_id_tot_ac[j] >= self.panel_B[0] and vol_id_tot_ac[j] <= self.panel_B[2]:
				AC_panel[j] = 'B'
				
				if vol_id_tot_ac[j] == self.panel_B[0]:
					AC_subpanel[j] = 1
					
				if vol_id_tot_ac[j] == self.panel_B[1]:
					AC_subpanel[j] = 2
				
				if vol_id_tot_ac[j] == self.panel_B[2]:
					AC_subpanel[j] = 3


			if vol_id_tot_ac[j] == self.panel_top:
				AC_panel[j] = 'T'
				AC_subpanel[j] = 0

	
		AC_panel = np.array(AC_panel)
		AC_subpanel = np.array(AC_subpanel)
		event_id_tot_ac = np.array(event_id_tot_ac)
		energy_dep_tot_ac = np.array(energy_dep_tot_ac)
		pair_flag_tot_ac = np.array(pair_flag_tot_ac)
		vol_id_tot_ac = np.array(vol_id_tot_ac)

		self.AC_panel = AC_panel 
		self.AC_subpanel = AC_subpanel 
		self.event_id_tot_ac = event_id_tot_ac 
		self.energy_dep_tot_ac = energy_dep_tot_ac
		self.pair_flag_tot_ac = pair_flag_tot_ac
		self.vol_id_tot_ac = vol_id_tot_ac
		
		
	######## COMPTON/PAIR AC EVENTS ######

	# IT SEPARATES THE COMPTON EVENTS FROM THE PAIR PRODUCTION ONES
	
    def compton_ac(self):

		where_compton_ac = np.where(self.pair_flag_tot_ac == 2)
		where_compton_ac = where_compton_ac[0]

		if len(where_compton_ac) != 0:
			event_id_tot_ac_compton = self.event_id_tot_ac[where_compton_ac]
			AC_panel_compton = self.AC_panel[where_compton_ac]
			AC_subpanel_compton = self.AC_subpanel[where_compton_ac]
			energy_dep_tot_ac_compton = self.energy_dep_tot_ac[where_compton_ac]
			pair_flag_tot_ac_compton = self.pair_flag_tot_ac[where_compton_ac]
  			
			self.event_id_tot_ac_compton = event_id_tot_ac_compton
			self.AC_panel_compton = AC_panel_compton 
			self.AC_subpanel_compton = AC_subpanel_compton 
			self.energy_dep_tot_ac_compton = energy_dep_tot_ac_compton 
			self.pair_flag_tot_ac_compton = pair_flag_tot_ac_compton 

	
		else:
			print('No Compton AC events found')
	
		self.where_compton_ac = where_compton_ac 

		
    def pair_ac(self):

		where_pair_ac = np.where(self.pair_flag_tot_ac == 1)
		where_pair_ac = where_pair_ac[0]
	

		if len(where_pair_ac) != 0:
			event_id_tot_ac_pair = self.event_id_tot_ac[where_pair_ac]
			AC_panel_pair = self.AC_panel[where_pair_ac]
			AC_subpanel_pair = self.AC_subpanel[where_pair_ac]
			energy_dep_tot_ac_pair = self.energy_dep_tot_ac[where_pair_ac]
			pair_flag_tot_ac_pair = self.pair_flag_tot_ac[where_pair_ac]
  			
			self.event_id_tot_ac_pair = event_id_tot_ac_pair
			self.AC_panel_pair = AC_panel_pair 
			self.AC_subpanel_pair = AC_subpanel_pair 
			self.energy_dep_tot_ac_pair = energy_dep_tot_ac_pair 
			self.pair_flag_tot_ac_pair = pair_flag_tot_ac_pair 


		else:
			print('No Pair AC events found')

		self.where_pair_ac = where_pair_ac

    def ray_ac(self):

		where_ray_ac = np.where(self.pair_flag_tot_ac == 3)
		where_ray_ac = where_ray_ac[0]
	

		if len(where_ray_ac) != 0:
			event_id_tot_ac_ray = self.event_id_tot_ac[where_ray_ac]
			AC_panel_ray = self.AC_panel[where_ray_ac]
			AC_subpanel_ray = self.AC_subpanel[where_ray_ac]
			energy_dep_tot_ac_ray = self.energy_dep_tot_ac[where_ray_ac]
			pair_flag_tot_ac_ray = self.pair_flag_tot_ac[where_ray_ac]
  			
			self.event_id_tot_ac_ray = event_id_tot_ac_ray
			self.AC_panel_ray = AC_panel_ray 
			self.AC_subpanel_ray = AC_subpanel_ray 
			self.energy_dep_tot_ac_ray = energy_dep_tot_ac_ray 
			self.pair_flag_tot_ac_ray = pair_flag_tot_ac_ray 


		else:
			print('No Rayleigh AC events found')

		self.where_ray_ac = where_ray_ac
