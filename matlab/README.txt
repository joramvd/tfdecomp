This repository contains both Matlab (working version) and Python (in progress) code for decomposing raw EEG data into time-frequency resolved data.

The backbone of this analysis is wavelet decomposition. FFT is applied on both the raw data, as well as wavelets of several frequencies, and multiplication is done in the frequency domain, after which the inverse FFT is taken.
This gives both power and phase estimates of frequency-specific activity over time.

With phase, the option exists to also do seeded connectivity analysis, or all-to-all connectivity analysis (e.g. for further graph-network analyses).
For connectivity, two metrics are implemented: inter-site phase clustering (or 'phase-locking value') and debiased weighted phase-lag index.

In the Matlab code, the raw EEG data needs to be in eeglab format. However, any format can be easily converted to eeglab format. For example, if you have Fieldtrip data resulting from ft_timelockanalysis or tf_preprocessing, there are simple conversion functions (fieldtrip2eeglab).

The Matlab package contains two scripts:
run_tfdecomp.m builds up a configuration structure that is passed on to tfdata = tfdecomp(cfg)
For basic power, phase and connectivity estimates for an eeglab structure with multiple experimental conditions, and several of those structures for a group of subjects, with stim-locked data (you can also relock the data to, e.g., the response using this function), there is probably no need to change the tfdecomp.m code. Just specify your ingredients in the run_tfdecomp.m script by filling in the cfg.
Type in help tfdecomp to see what is needed in the cfg.

For questions please e-mail me: joramvandriel@gmail.com