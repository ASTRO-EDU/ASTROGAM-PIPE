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
from DHSim import *

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
ene_min = sys.argv[8]              		 # Enter miminum energy [MeV]:
ene_max = sys.argv[9]               	 # Enter maximum energy [MeV]:
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

print astrogam_version

dhseASTROGAM = DHSim()

### parametri iniziali

dhseASTROGAM.parameter(sim_type, py_list, part_type, ene_range, pol_type, source_g, isStrip, repli, cal_flag, ac_flag, astrogam_version, bogemms_tag, ene_min, ene_max , passive_flag, energy_thresh)

#parte lettura file fits

dhseASTROGAM.get_filepath(astrogam_version, theta_type, N_in, part_type)

dhseASTROGAM.get_outdir(astrogam_version, theta_type, N_in, part_type, energy_thresh)


while ifile <= n_fits:
	start = time.time()
	
	print('Reading the THELSIM file.....'+ str(ifile))

	dhseASTROGAM.reading_fits(ifile, cal_flag, ac_flag, part_type, isStrip)
		
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#%                             Processing the tracker                          %
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if astrogam_version == 'V1.0' or astrogam_version == 'V1.1' or astrogam_version == 'V2.0' or astrogam_version == 'V10.0':

		# From Tracker volume ID to strip and tray ID and conversion from tray ID (starting from bottom) to plane ID (starting from the top)

		dhseASTROGAM.conversion(isStrip, repli)	

	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('                             Tracker   '                     )
	print('           Saving the Tracker raw hits (fits and .dat)      ')
	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	if astrogam_version == 'V1.0' or astrogam_version == 'V1.1' or astrogam_version == 'V2.0' or astrogam_version == 'V10.0':
		
		dhseASTROGAM.writing_G4raw(N_in, part_type, theta_type, phi_type, ifile, astrogam_version)		

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
			
			dhseASTROGAM.writing_AA_fake(astrogam_version, N_in, part_type, theta_type, phi_type, ifile)		
				
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('                             Tracker   ')
		print('  Creation of LUT data table with summed energy for each volume       ')
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

		# Loading the LUT

		if isStrip == 1:
			
			dhseASTROGAM.loading_geom(astrogam_version)
					
			### Flag events #####

			dhseASTROGAM.flag_events()

			###################
			# Accoppiamento Capacitivo
			###################
			
			dhseASTROGAM.acap()

			##################
			# Introducing the noise in every strip
			#################
			
			#dhseASTROGAM.noise()

			#
			# Summing the energy along the strip and applying the energy threshold
			#
				
			dhseASTROGAM.summing_energy(astrogam_version, N_in, part_type, theta_type, phi_type, ifile)
			
			#################
			# Energy Threshold
			################

			dhseASTROGAM.energy_threshold()		
		
			####index uniq event_id_tot

			dhseASTROGAM.index_uniq()

			
	if astrogam_version == 'V1.0' or astrogam_version == 'V1.1' or astrogam_version == 'V2.0' or astrogam_version == 'V10.0':

		if isStrip == 1:

			#Total number of strips

			print('Number of tracker triggered events: '+ str(dhseASTROGAM.N_trig))
			
			dhseASTROGAM.strip_analysis()
 		
			print('N_ev: '+ str(dhseASTROGAM.N_ev))				

			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('              Build the LEVEL 0 output            ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')		

			dhseASTROGAM.build_L0()

			# Level 0 = energy summed
			# Level 0 = the events are sorted in tray, and Y before X within the same tray
			# energy threshold applied		
			
			dhseASTROGAM.writing_L0(astrogam_version, N_in, part_type, theta_type, phi_type, ifile) 

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

			dhseASTROGAM.writing_AA_strip(N_in, part_type, ang_type, theta_type, phi_type, ifile)

			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('       L0.5 - cluster baricenter ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

			print('N_trig: '+ str(dhseASTROGAM.N_trig))

			
			dhseASTROGAM.L05_cluster_x()

			dhseASTROGAM.L05_cluster_y()

			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('                      Tracker   ')
			print('             L0.5 - X-Y layers merging ')
			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

			dhseASTROGAM.merging()			

			# Level 0.5 = energy summed, MIP threshold applied, strip position used

			dhseASTROGAM.writing_L05(astrogam_version, N_in, part_type, theta_type, phi_type, ifile)


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

			dhseASTROGAM.writing_AA_cluster(N_in, part_type, ang_type, theta_type, phi_type, ifile)
					
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

			dhseASTROGAM.writing_AA_cluster_pairs(N_in, part_type, ang_type, theta_type, phi_type, ifile)

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
									   
			dhseASTROGAM.writing_AA_cluster_compton(N_in, part_type, ang_type, theta_type, phi_type, ifile)

			print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			print('         ASCII data format for AA input - rayleigh')
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
			# - c11 = pair flag (3 = rayleigh)
									   
			dhseASTROGAM.writing_AA_cluster_rayleigh(N_in, part_type, ang_type, theta_type, phi_type, ifile)
			
			
			####### S1 data format ##############
			# ID Type Edep VolumeID X Y Z
			# - ID = Event ID
			# - Type = event flag
			# - Edep = energy deposited in the strip/calorimeter bar
			# - volume ID = unique ID for the volume
			# - X, Y, Z = position of the center of the volume

			dhseASTROGAM.writing_S1_trk(N_in, part_type, ang_type, theta_type, phi_type, ifile)


			####################
			#Angular resolution
			###################
			
			#dhseASTROGAM.resolution(N_in, part_type, ang_type, theta_type, phi_type, ifile)

	if cal_flag == 1:
		
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		#                             Processing the calorimeter
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('          Calorimeter Bar Energy attenuation'        )                
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('                   Calorimeter                ')
		print('              Applying the minimum cut                ')
		print('                Summing the energy                ')
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

		dhseASTROGAM.writing_cal_raw(N_in, part_type, theta_type, phi_type, ifile, astrogam_version)

		dhseASTROGAM.G4_cal()
		
		dhseASTROGAM.writing_G4cal(astrogam_version, N_in, part_type, theta_type, phi_type, ifile)

		dhseASTROGAM.compton_cal()
		
		if len(dhseASTROGAM.where_compton_cal) != 0:
			dhseASTROGAM.writing_G4_cal_compton(astrogam_version, N_in, part_type, theta_type, phi_type, ifile) 
		
		dhseASTROGAM.pair_cal()
		
		if len(dhseASTROGAM.where_pair_cal) != 0:
			dhseASTROGAM.writing_G4_cal_pair(astrogam_version, N_in, part_type, theta_type, phi_type, ifile) 

		dhseASTROGAM.ray_cal()
		
		if len(dhseASTROGAM.where_ray_cal) != 0:
			dhseASTROGAM.writing_G4_cal_ray(astrogam_version, N_in, part_type, theta_type, phi_type, ifile) 
		
		dhseASTROGAM.cal_sum()
		
		dhseASTROGAM.writing_cal_sum(astrogam_version, N_in, part_type, theta_type, phi_type, ifile)
		
		dhseASTROGAM.writing_S1_cal(N_in, part_type, ang_type, theta_type, phi_type, ifile)

		

	if ac_flag == 1:

		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('                          AC')
		print('                    write raw data                 ')
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

		dhseASTROGAM.writing_ac_raw(N_in, part_type, theta_type, phi_type, ifile, astrogam_version)

		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print('                          AC')
		print('                  Summing the energy                ')
		print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

		dhseASTROGAM.AC_analysis()

		dhseASTROGAM.writing_G4ac(astrogam_version, N_in, part_type, theta_type, phi_type, ifile)
		
		dhseASTROGAM.writing_S1_ac(N_in, part_type, ang_type, theta_type, phi_type, ifile)

		dhseASTROGAM.compton_ac()
		
		if len(dhseASTROGAM.where_compton_ac) != 0:
			dhseASTROGAM.writing_G4_ac_compton(astrogam_version, N_in, part_type, theta_type, phi_type, ifile) 
		
		dhseASTROGAM.pair_ac()
		
		if len(dhseASTROGAM.where_pair_ac) != 0:
			dhseASTROGAM.writing_G4_ac_pair(astrogam_version, N_in, part_type, theta_type, phi_type, ifile) 

		dhseASTROGAM.ray_ac()
		
		if len(dhseASTROGAM.where_ray_ac) != 0:
			dhseASTROGAM.writing_G4_ac_ray(astrogam_version, N_in, part_type, theta_type, phi_type, ifile) 
			



############    next fits file   ####################

	end = time.time() - start

	print('Time: ' + str(end))
				
	ifile = ifile + 1



	













