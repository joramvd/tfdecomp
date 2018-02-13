function [tf_pow, varargout] = tfdecomp(cfg,eegdat,varargin)

% Function for wavelet-based time-frequency decomposition of M/EEG data.
% Largely based on custom-written code from Mike X Cohen
%
% This function as it currently stands, needs a cell array with one cell
% per experimental condition, and each cell containing a channel-time-trial
% matrix of raw EEG data; this needs to be the second input argument.
%
% First input argument is a cfg structure (analogous to Fieldtrip functionality),
% with analysis-specific ingredients (see below).
% 
% Output arguments can be specified as follows (and largely depend on the
% cfg):
%
% - Power only:
% [tfpow,dim] = tfdecomp(cfg,eegdat);
% tfpow: 4D matrix of condition-channel-freq-time
% dim: structure containing time,freq,channel info for plotting
% 
% - Power and phase:
% [tfpow,tfphase,dim] = tfdecomp(cfg,eegdat);
% tfphase: similar to tfpow
%
% - Power and connectivity:
% [tfpow,tfphase,tfsync,dim] = tfdecomp(cfg,eegdat);
% tfsync: 5D[/6D] matrix of condition-seed-channel,freq,time[,metric], with
% metric being ISPC and/or dwPLI 
% 
% - Connectivity only:
% [tfsync,dim] = tfdecomp(cfg,eegdat);
% 
% - Power and regression weights:
% [tfpow,tfphase,tfbvals,dim] = tfdecomp(cfg,eegdat,regressors);
% tfbvals: 5D matrix of condition-channel-freq-time-regressors 
%
% - Regression weights and connectivity (no power):
% [tfbvals,tfsync,dim] = tfdecomp(cfg,eegdat,regressors);
% 
% The cfg should contain:
%
% -- Path/filenames for saving:
% cfg.writdir = 'path/to/folder/;
% cfg.filename = 'example_tf.mat';
%
% -- raw data specifics needed to compute filter ingredients:
% -- the sampling rate is either contained in the eeglab or fieldtrip data
% -- strucures, or you happen to know the specific number, just make sure it
% -- corresponds to the data used for input
% -- eegtime has all the time points of the epoch, in milliseconds
% -- (Fieldtrip uses seconds so in that case you need to converse)
% cfg.srate = EEG.srate; 
% cfg.eegtime = EEG.times; 
%
% -- if you relocked data to a response or other event, a stim-locked 
% -- baseline is still recommended, which requires saved time points of 
% -- when stimulus was presented at each trial:
% cfg.relock = baselinepoints; % see example run_tfdecomp.m on Github page,
% or email me: joramvandriel@gmail.com
%
% -- channel info: number of channels to analyze, and which channel label
% -- to take as seed for connectivity analysis
% cfg.channels = 1:64;
% cfg.connectivity = 'both'; % 'pli','iscp','both','none'
% cfg.seeds = {'fcz'}; % leave empty ({}) if no connectivity
%
% -- robust regression: fit a regression model with power as Y and e.g. RT
% -- as X, in Y=a+bX+e; the b-values will be stored for output; this requires
% -- the regressors (cell array with cell per condition, each cell with
% -- trial-by-regressors matrix) as a third input argument
% cfg.robfit = true; 
%
% -- frequency, time and baseline info:
% cfg.frequencies = [2 40 25]; % from min to max in nsteps
% cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
% cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled
% cfg.times2save = -200:25:1000;
% cfg.basetime = [-500 -200];
% cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
% cfg.erpsubract = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
% cfg.matchtrialn = true; % if true, conditions will be equated in terms of trial count (so SNR is comparable across conditions)
%
% -- other administrative stuff:
% cfg.report_progress = true;
% cfg.save_output = false;
% cfg.overwrite = false;

%% unpack cfg and setup settings

v2struct(cfg)

% if not specified, set fields to false
if isfield(cfg,'relock')
    startpoint = relock;
    relock = true;
else
    relock = false;
end

fields2check = {'basetime','singletrial','erpsubtract','matchtrialn','connectivity','report_progress','overwrite','robfit'};
for fieldi = 1:length(fields2check)
    if isfield(cfg,fields2check(fieldi)) 
        if strcmp(eval(fields2check{fieldi}),'none')
            eval([fields2check{fieldi} ' = false;']);
        end
    else
        eval([fields2check{fieldi} ' = false;']);
    end
end

% check output dir
if ~exist(writdir,'dir')
    mkdir(writdir)
end
if strcmp(writdir(end),filesep)==0
    writdir = [writdir filesep];
end
disp(writdir)

% frequencies
if strcmp(scale,'log')
    frex=logspace(log10(frequencies(1)),log10(frequencies(2)),frequencies(3));
elseif strcmp(scale,'lin')
    frex=linspace(frequencies(1),frequencies(2),frequencies(3));
end
nfreqs = frequencies(3);

% gaussian width and time
ntimepoints = size(eegdat{1},2);
if strcmp(scale,'log')
    s=logspace(log10(cycles(1)),log10(cycles(2)),nfreqs)./(2*pi.*frex);
elseif strcmp(scale,'lin')
    s=linspace(cycles(1),cycles(2),nfreqs)./(2*pi.*frex);
end
t=-ntimepoints/srate/2:1/srate:ntimepoints/srate/2-1/srate;

wavelets = zeros(nfreqs,length(t));
for fi=1:nfreqs
    wavelets(fi,:)=exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*s(fi)^2));
end

ntrials = cellfun(@(x) size(x,3),eegdat);
nconds = length(eegdat);

%% start

outputfilename = [writdir filename];
disp(outputfilename)

%% optional ERP subraction to compute non-phase locked power

if erpsubtract
    eegdat = cellfun(@(x) x-repmat(mean(x,3), [1 1 size(x,3)]),eegdat,'UniformOutput',false);
end

%% optional regressor extraction for robustfit

if robfit || strcmp(connectivity,'wispc')
    regressors = varargin{1};
    nreg = size(regressors{1},2);
end

%% optional matching of conditions in trial count

if matchtrialn
    [nmin,cmin] = min(ntrials);
    for condi=1:nconds
        if condi~=cmin
            trialsel = randperm(ntrials(condi));
            eegdat{condi} = eegdat{condi}(:,:,sort(trialsel(1:nmin)));
            if relock
                startpoint{condi} = startpoint{condi}(sort(trialsel(1:nmin)));
            end
            if robfit || strcmp(connectivity,'wispc')
                regressors{condi} = regressors{condi}(sort(trialsel(1:nmin)),:);
            end
            ntrials(condi)=nmin;
        end
    end
end

%% initialize output matrices

% setup time indexing
times2saveidx = zeros(size(times2save));
for ti=1:length(times2save)
    [~,times2saveidx(ti)]=min(abs(eegtime-times2save(ti)));
end

% baseline time indices
if sum(basetime)
    
    [~,basetimeidx(1)] = min(abs(eegtime-basetime(1)));
    [~,basetimeidx(2)] = min(abs(eegtime-basetime(2)));
    
    if relock % relocking option requires a 'startpoint' vector point to the point in time that initially was time 0
        [~,minbase] = min(abs(eegtime-(eegtime(1)+trialpad)));
        for condi=1:nconds
            nobasetrial = zeros(1,ntrials(condi));
            for ei=1:ntrials(condi)
                if (basetimeidx(1) - startpoint{condi}(ei) + 1)<minbase % first baseline point should not be before one sec.; earlier timepoints fall in the window that we padded to remove edge artifacts
                    nobasetrial(ei)=1; % mark trial index that have this issue
                end
            end
            if sum(nobasetrial)>0
                eegdat{condi}(:,:,nobasetrial==1)=[];
                ntrials(condi) = size(eegdat{condi},3);
                startpoint{condi}(nobasetrial==1)=[];
                if robfit
                    regressors{condi}(nobasetrial==1,:)=[];
                end
            end
        end
    end
    baselinedata = zeros(nconds,length(channels),nfreqs);
end

% empty output matrices
if singletrial
    tf_pow = cell(size(eegdat));
    tf_phase = NaN;
elseif robfit
    tf_rvals = zeros(nconds,length(channels),nfreqs,length(times2saveidx),size(regressors{1},2));
else
    [tf_pow, tf_phase] = deal(zeros(nconds,length(channels),nfreqs,length(times2saveidx)));
end
if  strcmp(connectivity,'both')
    tf_sync  = zeros(nconds,length(seeds),length(channels),nfreqs,length(times2saveidx),2);
elseif strcmp(connectivity,'pli') || strcmp(connectivity,'ispc')
    tf_sync  = zeros(nconds,length(seeds),length(channels),nfreqs,length(times2saveidx));
elseif strcmp(connectivity,'wispc')
    tf_sync  = zeros(nconds,length(seeds),length(channels),nfreqs,length(times2saveidx),nreg);
else
    tf_sync = NaN;
end

%% Now decompose
reverseStr='';

% loop around conditions
for condi=1:nconds
    
    if strcmp(connectivity,'wispc')
        regreshaped = permute(repmat(regressors{condi},[1 1 length(frex) length(times2saveidx)]),[2 3 4 1]);
    end
    
    Lconv = pow2(nextpow2( ntimepoints*ntrials(condi) + ntimepoints-1 ));
    rawconv = zeros(length(channels),nfreqs,length(times2save),ntrials(condi));
    
    % loop around channels
    for chani=channels
        
        if report_progress
            % display progress
            msg = sprintf('Decomposing channel %i/%i of experimental condition %i/%i...',  chani,length(channels),condi,nconds);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
        
        try
            EEGfft = fft(reshape(eegdat{condi}(chani,:,:),1,ntimepoints*ntrials(condi)),Lconv);
        catch me
            warning('More channels specified than present in data, rest is set to zero')
            continue
        end
        
        % loop around frequencies
        for fi=1:nfreqs
            
            % convolve and get analytic signal
            m = ifft(EEGfft.*fft(wavelets(fi,:),Lconv),Lconv);
            m = m(1:(ntimepoints*ntrials(condi) + ntimepoints-1));
            m = reshape(m(floor((ntimepoints-1)/2):end-1-ceil((ntimepoints-1)/2)),ntimepoints,ntrials(condi));
            
            % populate
            rawconv(chani,fi,:,:) = m(times2saveidx,:);
            
            % baseline power
            if sum(basetime)
                if relock
                    for ei=1:size(m,2)
                        basetimeshift = basetimeidx - startpoint{condi}(ei) + 1;
                        baselinedata(condi,chani,fi) = baselinedata(condi,chani,fi) + mean(abs(m(basetimeshift(1):basetimeshift(2),ei)).^2);
                    end
                    baselinedata(condi,chani,fi) = baselinedata(condi,chani,fi)./ntrials(condi); % average across trials
                else
                    baselinedata(condi,chani,fi) = mean(mean(abs(m(basetimeidx(1):basetimeidx(2),:)).^2,1),2); % incase of no re-locking this is no issue
                end
            end
            
            % robust regression
            if robfit
                for ti=1:length(times2saveidx)
                    
                    [tempR,stats]=robustfit(regressors{condi},log(abs(m(times2saveidx(ti),:)').^2));
                    tf_rvals(condi,chani,fi,ti,:) = tempR(2:end)./stats.se(2:end);
                    if ti==1
                        w = warning('query','last');
                        warning('off',w.identifier);
                    end
                end % end time-loop
            end

            
        end % end frequency loop
    end % end channel loop
    if report_progress
        fprintf('done.\n')
        reverseStr='';
    end
    
    % inter-site connectivity
    if connectivity
        for chanx=1:length(seeds)
            
            if report_progress
                % display progress
                msg = sprintf('Syncing seed %i/%i of experimental condition %i/%i...',  chanx,length(seeds),condi,nconds);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            
            for chany=channels
                
                % cross-spectral density
                csdm = squeeze(rawconv(strcmpi(seeds{chanx},{chanlocs.labels}),:,:,:) .* conj(rawconv(chany,:,:,:)));
                
                % ISPC
                if strcmp(connectivity,'ispc') || strcmp(connectivity,'both')
                    tmpsync = abs(mean(exp(1i*angle(csdm)),3)); % note: equivalent to ispc(fi,:) = abs(mean(exp(1i*(angle(sig1)-angle(sig2))),2));
                end
                
                % weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
                if strcmp(connectivity,'pli') || strcmp(connectivity,'both')
                    imagsum      = sum(imag(csdm),3);
                    imagsumW     = sum(abs(imag(csdm)),3);
                    debiasfactor = sum(imag(csdm).^2,3);
                    tmppli = (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor);
                end
                
                if strcmp(connectivity,'wispc')
                    
                    realispc = zeros(length(frex),length(times2save),nreg);
                    for regi = 1:nreg
                        realispc(:,:,regi) = abs(mean( squeeze(regreshaped(regi,:,:,:)) .* exp(1i*angle(csdm)),3));
                    end
                    
                    % permutation testing
                    if ~isfield(cfg,'nperm')
                        nperm=500;
                    end
                    fakeispc = zeros(nperm,length(frex),length(times2save),nreg);
                    
                    for permi=1:nperm
                        for regi=1:nreg
                            fakeispc(permi,:,:,regi) = abs(mean( squeeze(regreshaped(regi,:,:,randperm(ntrials(condi)))) .* exp(1i*angle(csdm)),3));
                        end
                    end
                end
                
                
                if strcmp(connectivity,'ispc')
                    tf_sync(condi,chanx,chany,:,:) = tmpsync;
                elseif strcmp(connectivity,'pli')
                    tf_sync(condi,chanx,chany,:,:) = tmppli;
                elseif strcmp(connectivity,'both')
                    tf_sync(condi,chanx,chany,:,:,1) = tmpsync;
                    tf_sync(condi,chanx,chany,:,:,2) = tmppli;
                elseif strcmp(connectivity,'wispc')
                    tf_sync(condi,chanx,chany,:,:,:) = (realispc - squeeze(mean(fakeispc))) ./ squeeze(std(fakeispc));
                end
                
            end % chany
        end % chanx
        if report_progress
            fprintf('done.\n')
            reverseStr='';
        end
    end
    
    % single trial power
    if singletrial
        tf_pow{condi} = abs(rawconv).^2;
    elseif sum(basetime)
        tf_pow(condi,:,:,:)   = mean(abs(rawconv).^2,4); % raw power; baseline correction is done after the condition loop
        tf_phase(condi,:,:,:) = abs(mean(exp(1i*angle(rawconv)),4));
    end
    
end % end condition loop

%% db convert: condition-specific baseline
if sum(basetime)
    if strcmp(baselinetype,'conspec')
        tf_pow = 10*log10( tf_pow ./ repmat(baselinedata,[ 1 1 1 length(times2save) ]) );
    elseif strcmp(baselinetype,'conavg')
        tf_pow = 10*log10( tf_pow ./ repmat(mean(baselinedata,1),[ nconds 1 1 length(times2save) ]) );
    end
end

% for plotting outside main function
dim = [];
dim.times = times2save;
dim.freqs = frex;
try
    dim.chans = chanlocs(channels);
catch me
    warning('No channel locations known; possibly you are analyzing component time series instead of scalp-channel time series?')
end
if connectivity
    dim.seeds = seeds;
end
dim.ntrials = ntrials;
dim.cfg_prev = cfg;

% specifiy output arguments
if connectivity & ~robfit
    if sum(basetime)
        varargout{1} = tf_phase;
        varargout{2} = tf_sync;
        varargout{3} = dim;
    else
        tf_pow = tf_sync;
        varargout{1} = dim;
    end
elseif connectivity & robfit
    if sum(basetime)
        varargout{1} = tf_phase;
        varargout{2} = tf_sync;
        varargout{3} = tf_rvals;
        varargout{4} = dim;
    else
        tf_pow = tf_rvals;
        varargout{1} = tf_sync;
        varargout{2} = dim;
    end
elseif robfit & ~connectivity
    if sum(basetime)
        varargout{1} = tf_phase;
        varargout{2} = tf_rvals;
        varargout{3} = dim;
    else
        tf_pow = tf_rvals;
        varargout{1} = dim;
    end
else
    varargout{1} = tf_phase;
    varargout{2} = dim;
end

%% save results
if save_output
    if connectivity
        if sum(basetime)
            save(outputfilename,'tf_pow','tf_phase','tf_sync','dim');
        else
            save(outputfilename,'tf_sync','dim');
        end
    elseif singletrial
        save(outputfilename,'tf_pow','dim');
    elseif robfit
        if sum(basetime)
            save(outputfilename,'tf_rvals','tf_pow','tf_phase','dim');
        else
            save(outputfilename,'tf_rvals','dim');
        end
    else
        save(outputfilename,'tf_pow','tf_phase','dim');
    end
end

%%
end
