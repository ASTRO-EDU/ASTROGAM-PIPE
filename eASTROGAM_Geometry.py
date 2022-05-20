"""
 eASTROGAM_Geometry.py  -  description
 ---------------------------------------------------------------------------------
 building the look-up table for the e-ASTROGAM simulation analysis
 ---------------------------------------------------------------------------------
 copyright            : (C) 2017 V. Fioretti, G. Giannella, S. Guidotti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python eASTROGAM_Geometry.py <version> 
 example:
 python eASTROGAM_Geometry.py V1.1
 ---------------------------------------------------------------------------------
 Parameters:
 - version: e-ASTROGAM geometry version
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2017/12/14: bugs fixing
"""


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import pickle
import sys
from astropy.table import Table

sys.argv[0] = 'eASTROGAM_Geometry.py'

astrogam_version = sys.argv[1]           # Enter eASTROGAM release (e.g. V1.0):


outdir = './conf/'
print('Configuration files path: '+ outdir)

if not os.path.exists('./conf/'):
	out_dir = os.makedirs(outdir,0777)

theta_deg_point = 30.
phi_deg_point = 225.
theta_point = theta_deg_point*(math.pi/180.)
phi_point = phi_deg_point*(math.pi/180.)

theta_deg_plane = 30.0
phi_deg_plane = 0.
theta_plane = theta_deg_plane*(math.pi/180.)
phi_plane = phi_deg_plane*(math.pi/180.)

theta_deg_pol = 90.
phi_deg_pol = 20.
theta_pol = theta_deg_pol*(math.pi/180.)
phi_pol = phi_deg_pol*(math.pi/180.)

#se theta_deg_pol = 90 inserire 0 al posto cos(theta_pol)
pol_vec = np.array([math.sin(theta_pol)*math.cos(phi_pol), math.sin(theta_pol)*math.sin(phi_pol), math.cos(theta_pol)]) 

# valore cos idl  -4.3711390e-08
# valore theta pol idl 1.5707964
# valore acos 1.57079637051
#val theta pol 1.57079632679

M = np.array([[math.cos(theta_point)*math.cos(phi_point), math.cos(theta_point)*math.sin(phi_point), -math.sin(theta_point)], [-math.sin(phi_point), math.cos(phi_point), 0.], [math.sin(theta_point)*math.cos(phi_point), math.sin(theta_point)*math.sin(phi_point), math.cos(theta_point)]])

M_minus = np.linalg.inv(M)

pol_vec_n = np.dot(pol_vec, M)

pol_vec_new = []

j = 0
while j<len(pol_vec_n):
	pol_vec_n_r =  np.around(pol_vec_n[j]*1000000.)/1000000.	

	pol_vec_new.append(pol_vec_n_r)

	j = j + 1
	
x_vec_new = pol_vec_new[0]
y_vec_new = pol_vec_new[1]
z_vec_new = pol_vec_new[2]

if y_vec_new >= 0. and x_vec_new >= 0.: 
	theta_pol_global = (180./math.pi)*math.acos(z_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2. + z_vec_new**2.))
	phi_pol_global = (180./math.pi)*math.asin(y_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2.))

if y_vec_new >= 0. and x_vec_new < 0.:
	theta_pol_global = (180./math.pi)*math.acos(z_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2. + z_vec_new**2.))
	phi_pol_global = 180. - (180./math.pi)*math.asin(y_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2.))

if y_vec_new < 0. and x_vec_new == 0.:
	theta_pol_global = (180./math.pi)*math.acos(z_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2. + z_vec_new**2.))
	phi_pol_global = 360 + (180./math.pi)*math.asin(y_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2.))

if y_vec_new < 0. and x_vec_new < 0.:
	theta_pol_global = (180./math.pi)*math.acos(z_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2. + z_vec_new**2.))
	phi_pol_global = 180. - (180./math.pi)*math.asin(y_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2.))

if y_vec_new < 0. and x_vec_new > 0.:
	theta_pol_global = (180./math.pi)*math.acos(z_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2. + z_vec_new**2.))
	phi_pol_global = 360 + (180./math.pi)*math.asin(y_vec_new/math.sqrt(x_vec_new**2. + y_vec_new**2.))

theta_pol_global = theta_pol_global*(math.pi/180.)
phi_pol_global = phi_pol_global*(math.pi/180.)

# source height
h_s = 150.  #cm

# Global Geometry (V1.1):
#N_tray = 56  
#N_layer = 1
#N_strip = 3840
#pitch = 0.240   #mm
#Tray_side = 921.6  #mm

# Global Geometry (V10.0):
N_tray = 60  
N_layer = 1
N_strip = 760
pitch = 0.500   #mm
Tray_side = 380.  #mm


# Tracker geometry [mm]
Si_t = 0.500
tracker_pitch = 10.
Al_t = tracker_pitch - Si_t

dist_tray = 0.   #mm

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('%                eASTROGAM'+astrogam_version+'                    %')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('% - Number of trays:' + str(N_tray)) 
print('% - Number of strips:'+ str(N_strip))
print('% - Pitch [mm]:'+ str(pitch))
print('% - Tray side [mm]:'+ str(Tray_side)) 
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('% Tracker thicknesses:')
print('% - Silicon thickness [mm]:'+ str(Si_t))
print('% - Al honeycomb thickness [mm]:'+str(Al_t))
print('% ----------------------------------------------')
print('% - Plane pitch [mm]:'+ str(tracker_pitch))
print('% - Trays distance [mm]:'+ str(dist_tray))
print('% ----------------------------------------------')

Lower_module_t = 0.

z_start = 0.

Central_module_t = Al_t
Upper_module_t = Si_t

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('% Tracker heights:')
print('% - Lower module height [mm]:'+ str(Lower_module_t))
print('% - Central module height [mm]:'+ str(Central_module_t))
print('% - Upper module height [mm]:'+ str(Upper_module_t))
print('% - Tray height [mm]:'+ str(Lower_module_t + Central_module_t + Upper_module_t))


TRK_t = Lower_module_t + Central_module_t + Upper_module_t
k = 1
while k < N_tray:
	TRK_t = TRK_t + Lower_module_t + Central_module_t + Upper_module_t + dist_tray
	
	k = k + 1

z_end = TRK_t + z_start

print('% - Tracker height [mm]:'+ str(TRK_t))
print('% - Tracker Z start [mm]:'+ str(z_start))
print('% - Tracker Z end [mm]:'+ str(z_end))


# Total number of strips

Total_vol_x = N_tray*N_strip
Total_vol_y = N_tray*N_strip
 

Glob_vol_id_x_top = np.zeros(Total_vol_x)
Glob_pos_x_top = np.zeros(Total_vol_x) 
Glob_z_x_top = np.zeros(Total_vol_x)  
Glob_moth_id_x_top = np.zeros(Total_vol_x)  
Glob_Strip_id_x_top = np.zeros(Total_vol_x)  
Glob_Si_id_x_top = np.zeros(Total_vol_x)  
Glob_tray_id_x_top = np.zeros(Total_vol_x)  
Glob_plane_id_x_top = np.zeros(Total_vol_x)  
Glob_energy_dep_x_top = np.zeros(Total_vol_x)  
Glob_pair_flag_x_top = np.zeros(Total_vol_x)

Glob_vol_id_y_top = np.zeros(Total_vol_y)  
Glob_pos_y_top = np.zeros(Total_vol_y) 
Glob_z_y_top = np.zeros(Total_vol_y) 
Glob_moth_id_y_top = np.zeros(Total_vol_y) 
Glob_Strip_id_y_top = np.zeros(Total_vol_y) 
Glob_Si_id_y_top = np.zeros(Total_vol_y) 
Glob_tray_id_y_top = np.zeros(Total_vol_y) 
Glob_plane_id_y_top = np.zeros(Total_vol_y) 
Glob_energy_dep_y_top = np.zeros(Total_vol_y) 
Glob_pair_flag_y_top = np.zeros(Total_vol_y) 

# all strips are readout

# ----> X layer

Tray_t = Lower_module_t + Central_module_t + Upper_module_t

t = 0

while t < N_tray:
	
	LowerModulePos_z = t*(dist_tray + Tray_t) + (Lower_module_t/2.)

	UpperModulePos_z = t*(dist_tray + Tray_t) + Lower_module_t + Central_module_t + (Upper_module_t/2.)      
	pos_z_x_top = UpperModulePos_z -(Upper_module_t/2.) + Si_t/2.
	pos_z_y_top = UpperModulePos_z -(Upper_module_t/2.) + Si_t/2.


	#print(pos_z_x_top)
	copyM = 1000000 + 1000000*t

#non funziona
	s = 0
	while s < N_strip:
		
		index = (t*N_strip) + s

		Glob_moth_id_x_top[index] = copyM + 90000
		Glob_tray_id_x_top[index] = t+1
		Glob_plane_id_x_top[index] = N_tray - t
		Glob_Si_id_x_top[index] = 0
		Glob_Strip_id_x_top[index] = s
		Glob_energy_dep_x_top[index] = 0.
		Glob_vol_id_x_top[index] = s
		Glob_pair_flag_x_top[index] = 0
		Strip_pos_x_top = -(Tray_side/2.0) + (pitch/2.) + (pitch*s)
		Glob_pos_x_top[index] = Strip_pos_x_top/10.    #cm
		Glob_z_x_top[index] = pos_z_x_top/10.

		Glob_moth_id_y_top[index] = copyM + 90000
		Glob_tray_id_y_top[index] = t+1
		Glob_plane_id_y_top[index] = N_tray - t
		Glob_Si_id_y_top[index] = 1
		Glob_Strip_id_y_top[index] = s
		Glob_energy_dep_y_top[index] = 0.
		Glob_vol_id_y_top[index] = s
		Glob_pair_flag_y_top[index] = 0
		Strip_pos_y_top = -(Tray_side/2.0) + (pitch/2.) + (pitch*s)
		Glob_pos_y_top[index] = Strip_pos_y_top/10.    #cm
		Glob_z_y_top[index] = pos_z_y_top/10.
		
		s = s + 1

	t = t + 1

		
col1 = fits.Column(name='VOLUME_ID', format='J', array=Glob_vol_id_x_top)	
col2 = fits.Column(name='MOTHER_ID', format='J', array=Glob_moth_id_x_top)
col3 = fits.Column(name='TRAY_ID', format='I', array=Glob_tray_id_x_top)
col4 = fits.Column(name='PLANE_ID', format='I', array=Glob_plane_id_x_top)
col5 = fits.Column(name='TRK_FLAG', format='I', array=Glob_Si_id_x_top)
col6 = fits.Column(name='STRIP_ID', format='J', array=Glob_Strip_id_x_top)
col7 = fits.Column(name='XPOS', format='F20.5', array=Glob_pos_x_top)
col8 = fits.Column(name='ZPOS', format='F20.5', array=Glob_z_x_top)
col9 = fits.Column(name='E_DEP', format='F20.5', array=Glob_energy_dep_x_top)
col10 = fits.Column(name='PAIR_FLAG', format='I', array=Glob_pair_flag_x_top)

cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10])
tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
if os.path.exists(outdir+'ARCH.XSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS'):
	os.remove(outdir+'ARCH.XSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS')
	tbhdu.writeto(outdir+'ARCH.XSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS')
else:
	tbhdu.writeto(outdir+'ARCH.XSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS')

fits.setval(outdir+'ARCH.XSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS', 'COMMENT', value='Creator = V. Fioretti', ext=1)
fits.setval(outdir+'ARCH.XSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS', 'COMMENT', value='eASTROGAM release = '+astrogam_version, ext=1)





col1 = fits.Column(name='VOLUME_ID', format='J', array=Glob_vol_id_y_top)	
col2 = fits.Column(name='MOTHER_ID', format='J', array=Glob_moth_id_y_top)
col3 = fits.Column(name='TRAY_ID', format='I', array=Glob_tray_id_y_top)
col4 = fits.Column(name='PLANE_ID', format='I', array=Glob_plane_id_y_top)
col5 = fits.Column(name='TRK_FLAG', format='I', array=Glob_Si_id_y_top)
col6 = fits.Column(name='STRIP_ID', format='J', array=Glob_Strip_id_y_top)
col7 = fits.Column(name='YPOS', format='F20.5', array=Glob_pos_y_top)
col8 = fits.Column(name='ZPOS', format='F20.5', array=Glob_z_y_top)
col9 = fits.Column(name='E_DEP', format='F20.5', array=Glob_energy_dep_y_top)
col10 = fits.Column(name='PAIR_FLAG', format='I', array=Glob_pair_flag_y_top)

cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10])
tbhdu = fits.BinTableHDU.from_columns(cols)			
		
		
if os.path.exists(outdir+'ARCH.YSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS'):
	os.remove(outdir+'ARCH.YSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS')
	tbhdu.writeto(outdir+'ARCH.YSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS')
else:
	tbhdu.writeto(outdir+'ARCH.YSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS')

fits.setval(outdir+'ARCH.YSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS', 'COMMENT', value='Creator = V. Fioretti', ext=1)
fits.setval(outdir+'ARCH.YSTRIP.TOP.eASTROGAM'+astrogam_version+'.TRACKER.FITS', 'COMMENT', value='eASTROGAM release = '+astrogam_version, ext=1)


print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('% Output FITS files with X and Y strip positions')
print('% - ARCH.XSTRIP.TOP.ASTROGAM'+astrogam_version+'.TRACKER.FITS')
print('% - ARCH.YSTRIP.TOP.ASTROGAM'+astrogam_version+'.TRACKER.FITS')

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('% GPS Set-up for the point source position:')
print('% - theta [deg.]:'+ str(theta_deg_point))
print('% - phi [deg.]:'+ str(phi_deg_point))
print('% - source height [cm]:'+ str(h_s))
print('% ----------------------------------------------')

# tracker height
h_t = z_end/10. # cm

# source height respect to tracker
h_r = h_s - h_t

# source distance from (0,0)
radius = h_r*math.tan(theta_point)
x_s = ((math.cos(phi_point))*radius)
y_s = ((math.sin(phi_point))*radius)


P_x = -math.sin(theta_point)*math.cos(phi_point)
P_y = -math.sin(theta_point)*math.sin(phi_point)
P_z = -math.cos(theta_point)

print('% Point source position:')
print('% - X [cm]:'+ str(x_s))
print('% - Y [cm]:'+ str(y_s))
print('% - Z [cm]:'+ str(h_s))
print('% Point Source direction:')
print('% - P_x:'+ str(P_x))
print('% - P_y:'+ str(P_y))
print('% - P_z:'+ str(P_z))
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


# Plane center distance from (0,0)
radius_plane = h_r*math.tan(theta_plane)

# Plane center
c_x_plane = ((math.cos(phi_plane))*radius_plane)
c_y_plane = ((math.sin(phi_plane))*radius_plane)
halfx_plane = Tray_side/20.
halfy_plane = Tray_side/20.

# Plane Momenta

P_x = -math.sin(theta_plane)*math.cos(phi_plane)
P_y = -math.sin(theta_plane)*math.sin(phi_plane)
P_z = -math.cos(theta_plane)

print('% Plane (square) center position:')
print('% - X [cm]:'+ str(c_x_plane))
print('% - Y [cm]:'+ str(c_y_plane))
print('% - Z [cm]:'+ str(h_s))
print('% Plane (square) side position:')
print('% - Half X [cm]:'+ str(halfx_plane))
print('% - Half Y [cm]:'+ str(halfy_plane))
print('% Plane Source direction:')
print('% - P_x:'+ str(P_x))
print('% - P_y:'+ str(P_y))
print('% - P_z:'+ str(P_z))
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

# Polarimetry Momenta

P_x_pol = math.sin(theta_pol_global)*math.cos(phi_pol_global)
P_y_pol = math.sin(theta_pol_global)*math.sin(phi_pol_global)
P_z_pol = math.cos(theta_pol_global)

print('% Polarimetry Momenta:')
print('% - P_x:'+ str(P_x_pol))
print('% - P_y:'+ str(P_y_pol))
print('% - P_z:'+ str(P_z_pol))
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')







