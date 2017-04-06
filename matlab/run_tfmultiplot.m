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

%% Example of condition-average power plot
%  also gives you line plot of condition-specific power

cfg = [];

cfg.metric = {'pow'};
cfg.chan = {'po8','po4','o2'};
cfg.freq = [8 14];
cfg.time = [800 1600];
cfg.scale = 'log';
cfg.connames = {'condition A','condition B'};
cfg.markevents = [0 500 1500];

tfmultiplot(cfg,tf_pow,dim);

%% Example of condition-comparison connectivity plot
%  for data with four conditions, we here compare the average of condition
%  1 and 2 with the average of condition 3 and 4, seeded from FCz

cfg = [];

cfg.metric = 'pli';
cfg.seed = {'fcz'};
cfg.chan = {'p1','p2','pz'};
cfg.freq = [8 14];
cfg.time = [2000 2500];
cfg.scale = 'log';
cfg.connames = {'condition A+B > condition C+D'};
cfg.concomp = [1 2; 3 4];
cfg.markevents = [0 500 1500];

tfmultiplot(cfg,tf_sync,dim);

%% Example of condition average power, comparing groups of channels
%  useful for e.g. lateralization analysis

cfg = [];

cfg.metric = {'pow'};
cfg.chan = {'po8','po4','o2'};
cfg.chandiff = {'po7','po3','o1'};
cfg.freq = [8 14];
cfg.time = [500 1200];
cfg.scale = 'log';
cfg.connames = {'condition A','condition B'};
cfg.markevents = [0 500 1500];

tfmultiplot(cfg,tf_pow,dim);


