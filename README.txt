This repository contains Matlab code for decomposing raw EEG data into time-frequency power and/or other metrics like inter-trial and inter-site phase clustering.

Input is a cfg structure with ingredients, and a eegdat cell array of n conditions, each with a channel-time-trial matrix.

The backbone of this analysis is wavelet decomposition. FFT is applied on both the raw data, as well as wavelets of several frequencies, and multiplication is done in the frequency domain, after which the inverse FFT is taken.
This gives both power and phase estimates of frequency-specific activity over time.

With phase, the option exists to also do seeded connectivity analysis, or all-to-all connectivity analysis (e.g. for further graph-network analyses).
For connectivity, two metrics are implemented: inter-site phase clustering (or 'phase-locking value') and debiased weighted phase-lag index.
There is also an option to subtract the ERP from the single trial raw data, so you get “induced” (so total minus “evoked”) power.

One new feature is robust regression (using Matlab’s robustfit), this requires a third input containing the regressors (X). Each time-frequency-channel raw power value will be Y in Y = a + bX1 + … + bXn + e, where the b-values for each time-frequency-channel point and each regressor are saved as regression weights. Regressors could be, for example, reaction time at each trial, or some continuous stimulus dimension. Robust regression does not require baseline correction.

You can also save single-trial power, e.g. if you want to export these data to more sophisticated machine-learning packages like sklearn (Python).

The Matlab package contains an example script that calls the tfdecomp function:
In run_tfdecomp.m you need to load the EEG data of one participant. Depending on your experiment and question, you divide the data up in conditions, and put them into separate cells of a variable (say, “eegdat”). Then set up a configuration structure that is passed on to tfdata = tfdecomp(cfg,eegdat).

This package contains four other functions:

- tfreject: sometimes single-trial power outliers can have a huge impact on the group-average results; these are usually caused by artifacts that were not detected during preprocessing. This function does a quick tf-decomposition and detects trials with outliers. I recommend to inspect the result of this detection in the raw data using EEGLAB.

- tfclustperm_tf: cluster-based permutation testing over time and frequency; a statistical test of a (difference) score against zero, correcting for multiple comparison. Channels are pre-selected.

- tfclustperm_chan: cluster-based permutation testing over space; time-frequency window or cluster is pre-selected. This function requires a neighbor structure offered by Fieldtrip, to go over candidate clusters in the data depending on the channel layout.

- tfmultiplot: plotting function at the group level (you need to average over participants beforehand), giving you in one glance: a time-frequency plot of a selection of channels (averaged over channels or subtracting channels), a topographical map of a time-frequency window, and a line plot of frequency-specific power. For a condition-contrast, the difference score is plotted in the time-frequency map, and the lines of the separate conditions are plotted in the line plot.


Future commits will include run_* example scripts of the above functions, with a 5-subject example dataset.

For questions please e-mail me: joramvandriel@gmail.com