% Example code of how to call tfdecomp function
% -- work in progress is to have an option to save single trial power for regression analysis


%%
cfg = [];

cfg.eeglab_path = 'Z:\Toolboxes\eeglab12_0_2_3b';
cfg.readdir = 'Z:\Stuff\Git\tfdecomp\matlab\';
cfg.writdir = cfg.readdir;
cfg.filename = '*data.mat';
cfg.projectname = 'sample';
cfg.relocking = 'none'; % 'none' or event code value; can be number if button response; or e.g. 'saccade' in case of simultaneous eye-tracking

cfg.seeds = {'oz'}; % channel name or empty {}
cfg.channels = 1:64;
cfg.connectivity = 'both'; % 'pli','iscp','both','none'
cfg.frequencies = [2 40 25]; % from min to max in nsteps
cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled
cfg.times2save = -200:25:1000;
cfg.basetime = [-500 -200];
cfg.stimbase = true; % in case of relocking the data, do you want the baseline to be pre-stimulus?
cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
cfg.erpsubtract = true; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
cfg.matchtrialn = true; % if true, conditions will be equated in terms of trial count (so SNR is comparable across conditions)

cfg.srate = 512;
cfg.epochtime = [-1 1.5];

cfg.report_progress = true;
cfg.save_output = true;
cfg.overwrite = false;
cfg.plot_output = []; % if left empty (so below removed or commented out), no plot will be generated
cfg.plot_output.chan = {'fcz','fc1','fc2'};
cfg.plot_output.freq = [3 8];
cfg.plot_output.time = [50 500];
cfg.plot_output.connames = {'conA','conB'};
cfg.plot_output.save = true;

[tf_pow, tf_phase, tf_sync, dim] = tfdecomp(cfg);

