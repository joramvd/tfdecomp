This repository contains Matlab code for decomposing raw EEG data into time-frequency resolved data.

Input is a cfg structure with ingredients, and a eegdat cell array of n conditions, each with a channel-time-trial matrix.

The backbone of this analysis is wavelet decomposition. FFT is applied on both the raw data, as well as wavelets of several frequencies, and multiplication is done in the frequency domain, after which the inverse FFT is taken.
This gives both power and phase estimates of frequency-specific activity over time.

With phase, the option exists to also do seeded connectivity analysis, or all-to-all connectivity analysis (e.g. for further graph-network analyses).
For connectivity, two metrics are implemented: inter-site phase clustering (or 'phase-locking value') and debiased weighted phase-lag index.
There is also an option to subtract the ERP from the single trial raw data, so you get “induced” (so total minus “evoked”) power.

One new feature is robust regression (using Matlab’s robustfit), this requires a third input containing the regressors (X). Each time-frequency-channel raw power value will be Y in Y = a + bX1 + … + bXn + e, where the b-values for each time-frequency-channel point and each regressor are saved as regression weights. Regressors could be, for example, reaction time at each trial, or some continuous stimulus dimension. Robust regression does not require baseline correction. 

The Matlab package contains two scripts:
In run_tfdecomp.m you need to set up a configuration structure that is passed on to tfdata = tfdecomp(cfg,eegdat)

You need to manually loop over participants
For basic power, phase and connectivity estimates for an eeglab structure with multiple experimental conditions, and several of those structures for a group of subjects, with stim-locked data (you can also relock the data to, e.g., the response using this function), there is probably no need to change the tfdecomp.m code. Just specify your ingredients in the run_tfdecomp.m script by filling in the cfg.
Type in help tfdecomp to see what is needed in the cfg.

tfdecomp also has the option to produce single-subject plots, for data quality checks.
The package also contains a separate plotting function tfmultiplot, which requires the output of tfdecomp as input, together with a cfg that contains the specification of what you want to plot. It will produce several subplots including a time-frequency map of a requested (group of) electrode(s), a topographical map of a requested time-frequency window, and line plots that visualize frequency-specific time courses for one or multiple experimental conditions.

Currently it only accepts single-subject (or an average over a group of subjects) data. Group-level analysis plotting, including cluster-based permutation testing, is work in progress.

For questions please e-mail me: joramvandriel@gmail.com