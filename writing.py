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


###### G4_raw.fits  ########

def writing_G4raw(event_id_tr, vol_id_tr, moth_id_tr, tray_id, plane_id, Strip_id_x, Strip_id_y, en_dep_tr, x_en_tr, y_en_tr, z_en_tr, x_ex_tr, y_ex_tr, z_ex_tr, child_id_tr, proc_id_tr, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile):

		
	col1 = fits.Column(name='EVT_ID', format='I', array=event_id_tr)	
	col2 = fits.Column(name='VOL_ID', format='I', array=vol_id_tr)
	col3 = fits.Column(name='MOTH_ID', format='J', array=moth_id_tr)
	col4 = fits.Column(name='TRAY_ID', format='I', array=tray_id)
	col5 = fits.Column(name='PLANE_ID', format='I', array=plane_id)
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
		



#### L0.fits  ############

def writing_L0(Glob_event_id_test, Glob_vol_id_test, Glob_moth_id_test, Glob_tray_id_test, Glob_plane_id_test, Glob_Si_id_test, Glob_Strip_id_test, Glob_pos_test, Glob_zpos_test, Glob_energy_dep_test, Glob_pair_flag_test, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile, N_trig):


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
	

############ L05.fits  #####################

def writing_L05(Glob_event_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_Si_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_pair_flag_cluster, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile, N_trig):

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
			
			
#############  G4_Cal.fits  #############

def writing_G4cal(event_id_tot_cal, bar_id_tot, bar_ene_tot, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile):

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

			
######### Cal Sum.fits  ##################

def writing_cal_sum(event_id_tot_cal_sum, bar_ene_tot_sum, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile):

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


########## G4_AC.fits #################

def writing_G4ac(event_id_tot_ac, AC_panel, AC_subpanel, energy_dep_tot_ac, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile):

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


####### AA Strip ##############

def writing_AA_strip(Glob_event_id_test, Glob_Si_id_test, Glob_tray_id_test, Glob_plane_id_test, Glob_Strip_id_test, Glob_pos_test, Glob_zpos_test, Glob_energy_dep_test, Glob_pair_flag_test, outdir, sim_tag, N_in, part_type, sname, ene_dis, ang_type, ene_type, theta_type, phi_type, pol_string, ifile):

	if os.path.exists(outdir+'/'+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
		os.remove(outdir+'/'+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
		data = open(outdir+'/'+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
	else:
		data = open(outdir+'/'+sim_tag+'_STRIP_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')


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


##### AA cluster #####

def writing_AA_cluster(Glob_event_id_cluster, Glob_Si_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_pair_flag_cluster, Glob_Strip_number_cluster, outdir, sim_tag, N_in, part_type, sname, ene_dis, ang_type, ene_type, theta_type, phi_type, pol_string, ifile):

	if os.path.exists(outdir+'/'+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
		os.remove(outdir+'/'+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
		data = open(outdir+'/'+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
	else:
		data = open(outdir+'/'+sim_tag+'_CLUSTER_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
		

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
	

##### AA cluster pairs ###############

def writing_AA_cluster_pairs(Glob_event_id_cluster, Glob_Si_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_pair_flag_cluster, Glob_Strip_number_cluster, outdir, sim_tag, N_in, part_type, sname, ene_dis, ang_type, ene_type, theta_type, phi_type, pol_string, ifile):

	if os.path.exists(outdir+'/'+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
		os.remove(outdir+'/'+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
		data = open(outdir+'/'+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
	else:
		data = open(outdir+'/'+sim_tag+'_CLUSTER_PAIR_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')

	
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



######## AA cluster compton ##############

def writing_AA_cluster_compton(Glob_event_id_cluster, Glob_Si_id_cluster, Glob_tray_id_cluster, Glob_plane_id_cluster, Glob_pos_cluster, Glob_zpos_cluster, Glob_energy_dep_cluster, Glob_pair_flag_cluster, Glob_Strip_number_cluster, outdir, sim_tag, N_in, part_type, sname, ene_dis, ang_type, ene_type, theta_type, phi_type, pol_string, ifile):

	if os.path.exists(outdir+'/'+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
		os.remove(outdir+'/'+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
		data = open(outdir+'/'+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
	else:
		data = open(outdir+'/'+sim_tag+'_CLUSTER_COMPTON_'+str(N_in)+part_type+'_'+sname+'_'+ene_dis+'_'+ang_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')


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
	
	
######### AA Fake #########

def writing_AA_fake(event_id_tr, plane_id, x_pos, y_pos, z_pos, en_dep_tr, child_id_tr, proc_id_tr, outdir, astrogam_version, py_name, sim_name, stripname, sname, N_in, part_type, ene_type, theta_type, phi_type, pol_string, ifile):

	if os.path.exists(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'.'+py_name+'.'+sim_name+'.'+stripname+'.'+sname+'.'+str(N_in)+part_type+'.'+ene_type+'MeV.'+str(theta_type)+'.'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat'):
		os.remove(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat')
		data = open(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
	else:
		data = open(outdir+'/AA_FAKE_eASTROGAM'+astrogam_version+'_'+py_name+'_'+sim_name+'_'+stripname+'_'+sname+'_'+str(N_in)+part_type+'_'+ene_type+'MeV_'+str(theta_type)+'_'+str(phi_type)+'.'+pol_string+str(ifile)+'.dat', 'w')
	
	j = 0
	while 1:
		where_event_eq = np.where(event_id_tr == event_id_tr[j])
		where_event_eq = where_event_eq[0]

		plane_id_temp = plane_id[where_event_eq]
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

