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
from IPython import embed as dbstop

class tfdecomp(object):

	def __init__(self, data = {}, settings = {}):
 
		self.data  = data
		self.freqs = settings['freqs']
		self.times = settings['times']
		self.srate = settings['srate']


	def setup_wavelets(self):

		if self.freqs['scale'][1]=='o':
			frex = np.logspace(np.log10(self.freqs['min']),np.log10(self.freqs['max']),self.freqs['num'])
		elif self.freqs['scale'][1]=='i':
			frex = np.linspace(self.freqs['min'],self.freqs['max'],self.freqs['num'])

		dbstop()

		# gaussian width and time
		s = np.logspace(np.log10(self.freqs['cycles'][0]),np.log10(self.freqs['cycles'][1]),self.freqs['num'])/(2*np.pi*frex)

		ntimepoints = self.data.shape[1]
		t = np.arange(-ntimepoints/self.srate/2., ntimepoints/self.srate/2., 1/self.srate)
		wavelets = np.empty([self.freqs['num'],ntimepoints])
		for f,fi in enumerate(frex):
        	wavelets[fi,:]=np.exp(2*1j*np.pi*fi*t).*np.exp(-t^2/(2*s[f]^2))











