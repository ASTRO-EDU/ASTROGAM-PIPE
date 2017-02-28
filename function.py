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



########## index_uniq da utilizzare solo con array gia' ordinati ########
def index_uniq(a):							#
									#
	a_index = []							#
									#
	for b in range(len(a)):						#
		a_index_old = b						#
		a_index.append(a_index_old)				#
									#
	a_uniq_index = []						#
									#
	b = 1								#
	while b < len(a):						#
									#
		if a[b] != a[b-1]:					#
			a_uniq_index.append(a_index[b-1])		#
									#
		if b == len(a)-1:					#
			a_uniq_index.append(a_index[b])			#
									#
		b = b + 1						#
									#
	return(a_uniq_index)						#
									#
#########################################################################






