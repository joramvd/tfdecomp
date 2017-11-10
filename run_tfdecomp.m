% Example code of how to call tfdecomp function

load sampleEEGdata.mat

% below is to illustrate how conditions are organized; it just takes the
% first and second half of all trials
eegdat{1} = double(EEG.data(:,:,1:floor(EEG.trials/2)));
eegdat{2} = double(EEG.data(:,:,floor(EEG.trials/2)+1:end));


%% for power, phase and phase-based connectivity

cfg = [];

cfg.seeds = {'fcz'}; % channel label; connectivity will be between this electrode and every other electrode (including the seed itself; for ispc this will be 1.0, for dwpli this will be NaN)
cfg.channels = 1:64;
cfg.chanlocs = EEG.chanlocs;
cfg.connectivity = 'pli'; % 'pli','ispc','both','none'
cfg.times2save = -500:25:1000;
cfg.frequencies = [1 40 25]; % from min to max in nsteps
cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled
cfg.basetime = [-500 -200]; % pre-stim baseline
cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
cfg.erpsubtract = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
cfg.matchtrialn = true; % if true, conditions will be equated in terms of trial count (so SNR is comparable across conditions)
cfg.srate = EEG.srate;
cfg.eegtime = EEG.times;

cfg.report_progress = true;
cfg.save_output = true;
cfg.overwrite = false;
cfg.writdir = pwd;
cfg.filename = 'tfdecomp_output.mat';

[tf_pow,tf_phase,tf_sync,dim] = tfdecomp(cfg,eegdat);


%% for all-to-all connectivity only

cfg = [];

cfg.seeds = {EEG.chanlocs(1:64).labels}; % give every channel label as seed
cfg.channels = 1:64;
cfg.chanlocs = EEG.chanlocs;
cfg.connectivity = 'ispc'; % 'pli','ispc','both','none'
cfg.times2save = -500:25:1000;
cfg.frequencies = [8 14 7]; % alpha band in 7 steps
cfg.cycles = [5 5]; % min max number of cycles used for min max frequency
cfg.scale = 'lin'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled

% omitting baseline specification will skip computing regular power; when
% another metric is specified (e.g. connectivity) then this will be the
% only output

cfg.erpsubtract = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
cfg.matchtrialn = true; % if true, conditions will be equated in terms of trial count (so SNR is comparable across conditions)
cfg.srate = EEG.srate;
cfg.eegtime = EEG.times;

cfg.report_progress = true;
cfg.save_output = true;
cfg.overwrite = false;
cfg.writdir = pwd;
cfg.filename = 'tfdecomp_output_graph.mat'; % all-to-all can be used for graph analysis

[tf_sync,dim] = tfdecomp(cfg,eegdat);

%% for robust regression

% robust regression: this requires a third input to the function: cell
% array with for each condition a regressor matrix of trial-by-nregressors
reg=zeros(EEG.trials,1);
for ei=1:EEG.trials
    [~,cue] = min(abs(cell2mat(EEG.epoch(ei).eventlatency)));
    reg(ei) = cell2mat(EEG.epoch(ei).eventlatency(cue+1));
end

regressors{1} = zscore(reg(1:floor(EEG.trials/2)));
regressors{2} = zscore(reg(floor(EEG.trials/2)+1:end));

cfg = [];

cfg.channels = 1:64;
cfg.chanlocs = EEG.chanlocs;
cfg.times2save = -500:25:1000;
cfg.frequencies = [2 80 50]; % 2 to 80 Hz in 50 steps
cfg.cycles = [5 5]; % min max number of cycles used for min max frequency
cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled

% omitting baseline specification will skip computing regular power; when
% another metric is specified (e.g. regression weights) then this will be the
% only output

cfg.erpsubtract = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
cfg.matchtrialn = true; % if true, conditions will be equated in terms of trial count (so SNR is comparable across conditions)
cfg.srate = EEG.srate;
cfg.eegtime = EEG.times;

cfg.report_progress = true;
cfg.save_output = true;
cfg.overwrite = false;
cfg.writdir = pwd;
cfg.filename = 'tfdecomp_output_robfit.mat'; % all-to-all can be used for graph analysis

cfg.robfit = true;

[tf_rvals,dim] = tfdecomp(cfg,eegdat,regressors);

%% for data re-locked to RT

% if you resorted your single trial data to be locked to another event than
% initially locked to during epoching, you should have another cell array
% for each condition, specifying for each trial what the new time point is
% of the stimulus (e.g. if RT-locking, the stimulus used to be at time 0
% but is now earlier in time, and this varies from trial to trial because
% now RT is always time 0 and the stim-RT interval differs from trial to
% trial)
% if you specify this, tfdecomp will compute single-trial baseline power
% relative to these baselinepoints, so you will have a pre-stimulus
% baseline
% one cautionary note: because the baseline shifts back in time, you need
% quite wide epochs because you still need zome buffer for edge artifacts
% for this reason you also need to specify the amount of time you padded to
% the trials during epoching to avoid edge artifacts; tfdecomp will check
% for each trial whether computig a baseline is still sensible, given this
% amount of trialpadding
% suppose we constructed epochs of -1.5 to 2 seconds, because we are
% interested in -.5 to 1 seconds, and we assumed a padding of 1 second is
% enough to avoid edge artifacts, the baseline should not fall inside this
% 1 second window, because according to our reasoning it will then be 
% contaminated by edge artifacts (especially for low frequencies); this 
% means that any RT of >500 ms will not give reliable results and should 
% thus be skipped.
% the example EEG data provided in this package is not ideal for this
% analysis because epochs start at -1 sec, but for illustration let's just
% say that the baseline can go all the way to -1 sec; this still means that
% trials with RT>500 ms cannot even be computed at all (because the
% baseline will start <1sec, which falls outside the epoch). This
% illustrates that during preprocessing, it's better to have wide epochs
% and remove redundant data later.

startpoint = zeros(EEG.trials,1);

for ei=1:EEG.trials
    
    [~,cue] = min(abs(cell2mat(EEG.epoch(ei).eventlatency)));
    
    % find RT and stim position in time
    [~,rtloc] = min(abs(EEG.times - EEG.epoch(ei).eventlatency{cue+1}));
    [~,zeroloc] = min(abs(EEG.times - 0 ));
    
    % find start point and replace with RT-aligned data
    startpoint(ei)    = rtloc-zeroloc+1;
    EEG.data(:,1:end-startpoint+1,ei) = EEG.data(:,startpoint:end,ei);
    
end
baselinepoints{1} = startpoint(1:floor(EEG.trials/2));
baselinepoints{2} = startpoint(floor(EEG.trials/2)+1:end);

cfg = [];

cfg.relock = baselinepoints;
cfg.trialpad = 0; % in ms

cfg.channels = 1:64;
cfg.chanlocs = EEG.chanlocs;
cfg.times2save = -500:25:1000;
cfg.frequencies = [1 40 25]; % from min to max in nsteps
cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled
cfg.basetime = [-500 -200]; % pre-stim baseline; NOTE: baseline remains stim-locked even when data are re-locked to RT (this is done with the baseline points saved above)
cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
cfg.erpsubtract = true; % NOTE: if true, phase-locked data in this case time-locked to the RT will be removed!
cfg.matchtrialn = true; % if true, conditions will be equated in terms of trial count (so SNR is comparable across conditions)
cfg.srate = EEG.srate;
cfg.eegtime = EEG.times;

cfg.report_progress = true;
cfg.save_output = true;
cfg.overwrite = false;
cfg.writdir = pwd;
cfg.filename = 'tfdecomp_output_rtlocked.mat';

[tf_pow,tf_phase,dim] = tfdecomp(cfg,eegdat); % take a look at the tf_phase output; with cfg.erpsubtract set to true, this shoul yield zero inter-trial phase clustering



