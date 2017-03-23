from tflevel1 import *

tf_settings = {}

# dirs
tf_settings['freqs'] = {}
tf_settings['freqs']['min'] = 1
tf_settings['freqs']['max'] = 40
tf_settings['freqs']['num'] = 25
tf_settings['freqs']['scale'] = 'log'
tf_settings['freqs']['cycles'] = [3, 12]

tf_settings['times'] = np.arange(-500,1500+25,25)
tf_settings['srate'] = 512

#data  = io.loadmat('sampleEEGdata.mat')
#eegdat = data['EEG'][0].data

data = np.random.rand(2,100,10) # two channels, 100 timepoints, 10 trials

check = tfdecomp(data, tf_settings)