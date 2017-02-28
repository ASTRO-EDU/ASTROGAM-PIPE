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




def conversion(vol_id_tr, isStrip, repli, moth_id_tr, tracker_bottom_vol_start, N_tray, tracker_top_bot_diff):

	# From Tracker volume ID to strip

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

	return(Strip_id_x, Strip_id_y, tray_id, plane_id)
