%% example code of how to call tfmultiplot

% input for tfmultiplot requires output of tfdecomp; this can be:
% - tf_pow (condition by channel by frequency by time matrix)
% - tf_phase (same dimensions as tf_pow)
% - tf_sync (additional dimension of seed electrode(s), and containing one
%   or two connectivity metrics [pli and/or ispc]
% - dim (structure containing timepoints, frequencies, and channel
%        locations)
%
% Currently this function only works for single-subject or group-averaged
% data
% You should carefully check if the dimensions of the input data match
% the requested metic to be plotted; look out for singleton-dimensions,
% as there is some squeezing done in the function.

%% Example of condition-average power plot
%  also gives you line plot of condition-specific power

cfg = [];

cfg.metric = {'pow'};
cfg.chan = {'po7','po3','o1'};
cfg.freq = [6 14];
cfg.time = [300 800];
cfg.scale = 'log';
cfg.connames = {'condition A'};
cfg.concomp = [2; 1];
cfg.markevents = [0 500];

tfmultiplot(cfg,tf_pow,dim);

%% Example of condition-comparison connectivity plot
%  for data with four conditions, we here compare the average of condition
%  1 and 2 with the average of condition 3 and 4, seeded from FCz

cfg = [];

cfg.metric = 'pli';
cfg.seed = {'fcz','fc1','fc2','fz'};
cfg.chan = {'p1','p2','pz','poz','po7','po8'};
cfg.freq = [2 7];
cfg.time = [50 400];
cfg.scale = 'log';
cfg.connames = {'condition A'};
% cfg.concomp = [1 2; 3 4];
cfg.markevents = [0];

tfmultiplot(cfg,squeeze(tf_sync(:,:,:,:,:,1)),dim);

%% Example of condition average power, comparing groups of channels
%  useful for e.g. lateralization analysis

cfg = [];

cfg.metric = {'pow'};
cfg.chan = {'po8','po4','o2'};
cfg.chandiff = {'po7','po3','o1'};
cfg.freq = [6 14];
cfg.time = [300 800];
cfg.scale = 'log';
cfg.connames = {'condition A'};
cfg.markevents = [0 500];

tfmultiplot(cfg,plotdat(1,:,:,:),dim);


