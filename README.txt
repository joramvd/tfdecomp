This repository contains both Matlab (working version) and Python (in progress) code for decomposing raw EEG data into time-frequency resolved data.

The backbone of this analysis is wavelet decomposition. FFT is applied on both the raw data, as well as wavelets of several frequencies, and multiplication is done in the frequency domain, after which the inverse FFT is taken.
This gives both power and phase estimates of frequency-specific activity over time.

With phase, the option exists to also do seeded connectivity analysis, or all-to-all connectivity analysis (e.g. for further graph-network analyses).
For connectivity, two metrics are implemented: inter-site phase clustering (or 'phase-locking value') and debiased weighted phase-lag index.

In the Matlab code, the raw EEG data needs to be in eeglab format. However, any format can be easily converted to eeglab format.

For questions please e-mail me: joramvandriel@gmail.com
