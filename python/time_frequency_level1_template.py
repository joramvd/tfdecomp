#%%
# Template level-1 time-frequency analysis script for Python
# Based on Matlab level-1 script by mikexcohen@gmail.com
# Joram van Driel, VU, January 2015

#%% Import some modules
import math
import glob
import scipy as sp
import numpy as np
from os import listdir
from os import chdir
from scipy import stats, io


class tfdecomp(object):

	def __init__(self, data = {}, settings = {}):
 
		self.data  = data
		self.freqs = settings['freqs']
		self.times = settings['times']
		self.freqscale = settings['freqscale']

	def data_relock(self):



	def setup_wavelets(self):

		if frequencyscaling[1]=='o':
			frex = np.logspace(np.log10(self.freqs['min']),np.log10(self.freqs['max']),self.freqs['num'])
		elif frequencyscaling[1]=='i':
			frex = np.linspace(self.freqs['min'],self.freqs['max'],self.freqs['num'])













