% Example code of how to call tfclustperm[_tf][_chan] functions

% Below is a heavily downsampled group data set of 5 participants, so the 
% file size is not too large
% It contains two conditions: error trials versus correct trials, and the
% data are time-locked to the button press (i.e. response-locked)

load sampleEEGdata.mat

%% Set up cfg

cfg = [];

cfg.type = 'T'; % at pixel-level: t-test or z-transform?

% This is the threshold for both the pixel- and cluster-level test.
% If not specified otherwise, it is used for a two-sided test (.025
% on either tail). If you wish a one-sided test, specify below but
% don't change the p-value (e.g. 0.05 will be used on one side only).
% At the cluster-level, the 95% CI threshold (1-p-val) is checked
% fot the maximum cluster-size under the null-hypothesis
% The pixel- and cluster-level thresholding are related in the following
% way: suppose you 
cfg.pval = 0.05; 
