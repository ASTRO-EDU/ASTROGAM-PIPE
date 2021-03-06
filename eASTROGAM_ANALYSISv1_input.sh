
#!/usr/bin/python

##############################
#     parametri iniziali     #
##############################

#astrogam_version			Enter eASTROGAM release (e.g. V1.0):
#bogemms_tag				Enter BoGEMMS release (e.g. 211):
#sim_type				Enter simulation type [0 = Mono, 1 = Range, 2 = Chen, 3: Vela, 4: Crab, 5: G400]:
#py_list				Enter the Physics List [0 = QGSP_BERT_EMV, 100 = ARGO, 300 = FERMI, 400 = ASTROMEV]:
#N_in					Enter the number of emitted particles:
#part_type				Enter the particle type [ph = photons, mu = muons, g = geantino, p = proton, el = electron]:
#ene_range				Enter energy distribution [0 = MONO, 1 = POW, 2 = EXP, 3 = LIN]:
#ene_min				Enter miminum energy [MeV]:
#ene_max				Enter maximum energy [MeV]:
#ang_type				Enter the angular distribution [e.g. UNI, ISO]:
#theta_type				Enter theta:
#phi_type				Enter phi:
#pol_type				Is the source polarized? [0 = false, 1 = true]:
#pol_angle				Enter Polarization angle:
#source_g				Enter source geometry [0 = Point, 1 = Plane]:
#isStrip				Strip/Pixels activated? [0 = false, 1 = true]
#repli					Strips/Pixels replicated? [0 = false, 1 = true]
#cal_flag				Is Cal present? [0 = false, 1 = true]:
#ac_flag				Is AC present? [0 = false, 1 = true]:
#passive_flag				Is Passive present? [0 = false, 1 = true]:
#energy_thresh				Enter energy threshold [keV]:
#ifile					Enter the initial number of FITS files
#n_fits					Enter the final number of FITS files:



#module load python2.7-sci


####per analizzare i file fits di boogems

python eASTROGAM_ANALYSISv1_file_remote.py V1.1 211 0 400 100000 ph 0 100 100 UNI 30 225 0 20 0 1 1 1 0 0 15 0 3


####per unire gli output processati da eASTROGAM_ANALYSISv3_file_remote

#python eASTROGAM_ANALYSISv1_all_remote.py V1.0 211 0 400 100000 ph 0 100 100 UNI 30 225 0 20 0 1 1 1 0 0 15 0 1

### per creare i LUT (file Geometry) ####

#python eASTROGAM_Geometry.py V1.1
