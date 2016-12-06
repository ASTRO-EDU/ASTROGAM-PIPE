from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math

#parametri iniziali

astrogam_version = 'V1.0'   # Enter eASTROGAM release (e.g. V1.0):
bogemms_tag = 211         # Enter BoGEMMS release (e.g. 211):
sim_type = 0            # Enter simulation type [0 = Mono, 1 = Range, 2 = Chen, 3: Vela, 4: Crab, 4: G400]:
py_list = 400          # Enter the Physics List [0 = QGSP_BERT_EMV, 100 = ARGO, 300 = FERMI, 400 = ASTROMEV]:
N_in = 100000                # Enter the number of emitted particles:
part_type = "ph"           # Enter the particle type [ph = photons, mu = muons, g = geantino, p = proton, el = electron]:
n_fits = 3              # Enter number of FITS files:
ene_range = 0          # Enter energy distribution [0 = MONO, 1 = POW, 2 = EXP, 3 = LIN]:
ene_min = .5             # Enter miminum energy [MeV]:
ene_max = 100            # Enter maximum energy [MeV]:
ang_type = "UNI"           # Enter the angular distribution [e.g. UNI, ISO]:
theta_type = 30         # Enter theta:
phi_type = 225           # Enter phi:
pol_type = 0           # Is the source polarized? [0 = false, 1 = true]: 
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

#sim_tag = 'eAST'+bogemms_tag+astrogam_tag+'0102'

if ene_range == 0:
	ene_dis = 'MONO'
	ene_type = ene_min
	if ene_type >= 1:
		ene_type = repr(ene_type)
	if ene_type < 1:
		ene_type = repr(ene_type)
		ene_type = ene_type[:]
	if type(ene_type) is int:
		print('It''s int')
	else:	
		nstring = len(ene_type)
		ene_type_notzero = ene_type
      		flag = 1 
		
		#for ichar_reverse in range(0, nstring):
			#ichar = (nstring-1) - ichar_reverse
		if ene_type[0] == '0' or  ene_type[0] == '.':
			if flag == 1:
				ene_type_notzero = ene_type_notzero[0]
		else:
			flag = 0
		ene_type = ene_type_notzero

print(ene_type)
	
	

















