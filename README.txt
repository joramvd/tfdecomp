This repository contains Matlab code for decomposing raw EEG data into time-frequency resolved data.

The backbone of this analysis is wavelet decomposition. FFT is applied on both the raw data, as well as wavelets of several frequencies, and multiplication is done in the frequency domain, after which the inverse FFT is taken.
This gives both power and phase estimates of frequency-specific activity over time.

With phase, the option exists to also do seeded connectivity analysis, or all-to-all connectivity analysis (e.g. for further graph-network analyses).
For connectivity, two metrics are implemented: inter-site phase clustering (or 'phase-locking value') and debiased weighted phase-lag index.
There is also an option to subtract the ERP from the single trial raw data, so you get “induced” (so total minus “evoked”) power.

In the Matlab code, the raw EEG data needs to be in eeglab format. However, any format can be easily converted to eeglab format. For example, if you have Fieldtrip data resulting from ft_timelockanalysis or tf_preprocessing, there are simple conversion functions (fieldtrip2eeglab).

The Matlab package contains two scripts:
In run_tfdecomp.m you need to set up a configuration structure that is passed on to tfdata = tfdecomp(cfg)
For basic power, phase and connectivity estimates for an eeglab structure with multiple experimental conditions, and several of those structures for a group of subjects, with stim-locked data (you can also relock the data to, e.g., the response using this function), there is probably no need to change the tfdecomp.m code. Just specify your ingredients in the run_tfdecomp.m script by filling in the cfg.
Type in help tfdecomp to see what is needed in the cfg.

tfdecomp also has the option to produce single-subject plots, for data quality checks.
The package also contains a separate plotting function tfmultiplot, which requires the output of tfdecomp as input, together with a cfg that contains the specification of what you want to plot. It will produce several subplots including a time-frequency map of a requested (group of) electrode(s), a topographical map of a requested time-frequency window, and line plots that visualize frequency-specific time courses for one or multiple experimental conditions.

Currently it only accepts single-subject (or an average over a group of subjects) data. Group-level analysis plotting, including cluster-based permutation testing, is work in progress.

For questions please e-mail me: joramvandriel@gmail.com