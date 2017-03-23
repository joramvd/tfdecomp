function [tf_pow, tf_phase, tf_sync, frex] = tfdecomp(cfg)

% Function for wavelet-based time-frequency decomposition of M/EEG data.
% Largely based on custom-written code from Mike X Cohen
%
% This function as it currently stands, needs a data structure according to the eeglab format
% It further needs a cfg structure (inspired by Fieldtrip functionality)
%
% This cfg necessitates the following:
%
% -- Pathnames to eeglab, data, filenames etc:
% cfg.eeglab_path = 'Z:\Toolboxes\eeglab12_0_2_3b';
% cfg.readdir = 'Z:\Stuff\Git\tfdecomp\matlab\';
% cfg.writdir = cfg.readdir;
% cfg.filename = '*data.mat';
%
% -- project name; this word, in addition to the first four letters of the
% -- filename (assuming a 'pp0x' prefix of 'proefpersoon number x') will be
% -- added to the outputfilename
% cfg.projectname = 'sample'; 
%
% -- raw data specifics needed to compute filter ingredients:
% cfg.srate = 256;
% cfg.epochtime = [-1 1.5];
% cfg.npoints = cfg.srate*sum(abs(cfg.epochtime));
%
% -- you can relock data to a response or other event:
% cfg.relocking = 'none'; % 'none' or event code value; can be number if button response; or e.g. 'saccade' in case of simultaneous eye-tracking 
%
% -- channel info: number of channels to analyze, and which channel label
% -- to take as seed for connectivity analysis
% cfg.channels = 1:64;
% cfg.connectivity = 'both'; % 'pli','iscp','both','none'
% cfg.seeds = {'fcz'}; % leave empty if no connectivity
%
% -- frequency, time and baseline info:
% cfg.frequencies = [2 40 25]; % from min to max in nsteps
% cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
% cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled
% cfg.times2save = -200:25:1000;
% cfg.basetime = [-500 -200];
% cfg.stimbase = true; % in case of relocking the data, do you want the baseline to be pre-stimulus?
% cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
% 
% -- number of conditions
% cfg.nconds = 1;
% 
% -- other administrative stuff:
% cfg.report_progress = true;
% cfg.save_output = false;
% cfg.plot_output.chan = {'poz','oz'};
% cfg.plot_output.freq = [8 16];
% cfg.plot_output.time = [200 700];
% cfg.plot_output.connames = {'con1'};
% cfg.save_plot = false;


%% unpack cfg and setup settings

v2struct(cfg)

path_list_cell = regexp(path,pathsep,'Split');
if (exist(eeglab_path,'dir') == 7) && ~any(ismember(eeglab_path,path_list_cell))
     % addpath(genpath(eeglab_path),'-begin');
     addpath(genpath(eeglab_path));
end

% files to use in analysis
filz = dir([ readdir filename ]);

% frequencies
if strcmp(scale,'log')
    frex=logspace(log10(frequencies(1)),log10(frequencies(2)),frequencies(3));
elseif strcmp(scale,'lin')
    frex=linspace(frequencies(1),frequencies(2),frequencies(3));
end
num_freqs = frequencies(3);

% gaussian width and time
s=logspace(log10(cycles(1)),log10(cycles(2)),num_freqs)./(2*pi.*frex);
t=-npoints/srate/2:1/srate:npoints/srate/2-1/srate;

wavelets = zeros(num_freqs,length(t));
for fi=1:num_freqs
    wavelets(fi,:)=exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*s(fi)^2));
end

%% loop through files

for subno=1:length(filz)
    
    outputfilename = [ writdir filz(subno).name(1:4) '_' projectname '_tfdecomp.mat' ];
    
    if exist(outputfilename,'file'), continue; end
    load([ readdir filz(subno).name ])
    
    % match ALLEEG to EEG if only one dataset
    if ~exist('ALLEEG','var'), ALLEEG=EEG; end;
    if ~exist('EEG','var'), EEG=ALLEEG(1); end;
    
    %% optional relocking of data
    %  e.g. if data were originally epoched around stim-onset, and you want
    %  to do response-locked analysis
    %  note: in case of e.g. slow responses and relocking to response, you
    %  need to make sure the epochs are wide enough, because the data are
    %  shifted but epoch length remains the same
    
    if ~strcmp(relocking,'none')
        
        % startpoint is a cell of vectors containing the real trial onset times (needed for baselining)
        startpoint    = cell(length(ALLEEG),1);
        baselinepoint = cell(length(ALLEEG),1);
        
        for condi=1:length(ALLEEG)
            
            startpoint{condi}    = zeros(length(ALLEEG(condi).epoch),1);
            baselinepoint{condi} = zeros(length(ALLEEG(condi).epoch),1);
            removetrial          = zeros(length(ALLEEG(condi).epoch),1);
            
            for ei=1:length(ALLEEG(condi).epoch)
                
                try
                    % find position
                    [~,reloc  ] = min(abs(ALLEEG(condi).times - ALLEEG(condi).epoch(ei).eventlatency{find(cell2mat(ALLEEG(condi).epoch(ei).eventlatency)==0)+1}));
                    [~,zeroloc] = min(abs(ALLEEG(condi).times - 0 ));
                    
                    % find start point and replace with re-aligned data
                    startpoint{condi}(ei)    = reloc-zeroloc+1;
                    baselinepoint{condi}(ei) = zeroloc - (reloc-zeroloc);
                    ALLEEG(condi).data(:,1:end-startpoint{condi}(ei)+1,ei) = ALLEEG(condi).data(:,startpoint{condi}(ei):end,ei);
                catch me, removetrial(ei)=1;
                end
            end
            
            % now remove trials with no response
            if sum(removetrial)>0
                fprintf('Note: %i additional trials are removed that could not be relocked to event of interest!\n', sum(removetrial));
                ALLEEG(condi) = pop_select(ALLEEG(condi),'notrial',find(removetrial));
            end
        end % end condition loop
    end
    
    %% initialize output matrices
    
    % setup time indexing
    times2saveidx = zeros(size(times2save));
    for ti=1:length(times2save)
        [~,times2saveidx(ti)]=min(abs(EEG.times-times2save(ti)));
    end
    
    % baseline time indices
    [~,basetimeidx(1)] = min(abs(EEG.times-basetime(1)));
    [~,basetimeidx(2)] = min(abs(EEG.times-basetime(2)));
    
    % empty output matrices
    [tf_pow, tf_phase] = deal(zeros(length(ALLEEG),length(channels),num_freqs,length(times2saveidx)));
    if  strcmp(cfg.connectivity,'both')
        tf_sync  = zeros(length(ALLEEG),length(seeds),length(channels),num_freqs,length(times2saveidx),2);
    elseif strcmp(cfg.connectivity,'pli') || strcmp(cfg.connectivity,'iscp')
        tf_sync  = zeros(length(ALLEEG),length(seeds),length(channels),num_freqs,length(times2saveidx));
    end
    
    baselinedata = zeros(length(ALLEEG),ALLEEG(1).nbchan,num_freqs);
    rawconv = cell(size(ALLEEG));
    
    %% Now decompose
    reverseStr='';

    % loop around conditions
    for condi=1:nconds
        
        Lconv = pow2(nextpow2( npoints*ALLEEG(condi).trials + npoints-1 ));
        
        % loop around channels
        for chani=channels
            
            if report_progress
                % display progress
                msg = sprintf('Decomposing channel %i/%i of experimental condition %i/%i...',  chani,length(channels),condi,nconds);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            
            EEGfft = fft(reshape(ALLEEG(condi).data(chani,:,:),1,npoints*ALLEEG(condi).trials),Lconv(condi));
            
            % loop around frequencies
            for fi=1:num_freqs
                
                % convolve and get analytic signal
                m = ifft(EEGfft.*fft(wavelets(fi,:),Lconv),Lconv);
                m = m(1:(npoints*ALLEEG(condi).trials + npoints-1));
                m = reshape(m(floor((npoints-1)/2):end-1-ceil((npoints-1)/2)),npoints,ALLEEG(condi).trials);

                % enter into TF matrix
                tf_pow(condi,chani,fi,:)   = mean(abs(m(times2saveidx,:)).^2,2)';
                tf_phase(condi,chani,fi,:) = abs(mean(exp(1i*angle(m(times2saveidx,:))),2))';
                
                % also need single-trial data for baselining
                if ~strcmp(relocking,'none')
                    for ti=1:size(m,2)
                        baselinedata(condi,chani,fi) = baselinedata(condi,chani,fi) + mean(abs(m(baselinepoint{condi}(ti)-min(baselinepoint{condi}(ti)-1,basetimeidx(1)):baselinepoint{condi}(ti)-basetimeidx(2),ti).^2));
                    end
                    baselinedata(condi,chani,fi) = baselinedata(condi,chani,fi)./ALLEEG(condi).trials; % average across trials
                else
                    baselinedata(condi,chani,fi) = mean(mean(abs(m(basetimeidx(1):basetimeidx(2),:)).^2,1),2);
                end
                
                % for inter-site connectivity we need raw convolution for cross-spectral density matrix
                if ~isempty(seeds)
                    % initialize
                    if chani==1 && fi==1
                        rawconv = zeros(length(channels),num_freqs,length(times2save),ALLEEG(condi).trials);
                    end
                    % populate
                    rawconv(chani,fi,:,:) = m(times2saveidx,:);
                end
                
            end % end frequency loop
        end % end channel loop
        if report_progress
            fprintf('done.\n')
            reverseStr='';
        end

        % inter-site connectivity
        if ~isempty(seeds)
            for chanx=1:length(seeds)
                
                if report_progress
                    % display progress
                    msg = sprintf('Syncing seed %i/%i of experimental condition %i/%i...',  chanx,length(seeds),condi,nconds);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
                
                for chany=channels
                    
                    % cross-spectral density
                    csd = squeeze(rawconv(strcmpi(seeds{chanx},{EEG.chanlocs.labels}),:,:,:) .* conj(rawconv(chany,:,:,:)));
                    
                    % ICPS
                    if strcmp(cfg.connectivity,'icps') || strcmp(cfg.connectivity,'both')
                        tmpsync = abs(mean(exp(1i*angle(csd)),3)); % note: equivalent to ispc(fi,:) = abs(mean(exp(1i*(angle(sig1)-angle(sig2))),2));
                    end
                    
                    % weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
                    if strcmp(connectivity,'pli') || strcmp(connectivity,'both')
                        imagsum      = sum(imag(csd),3);
                        imagsumW     = sum(abs(imag(csd)),3);
                        debiasfactor = sum(imag(csd).^2,3);
                        tmppli = (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor);
                    end
                    
                    if strcmp(connectivity,'icps')
                        tf_sync(condi,chanx,chany,:,:) = tmpsync;
                    elseif strcmp(connectivity,'pli')
                        tf_sync(condi,chanx,chany,:,:) = tmppli;
                    elseif strcmp(connectivity,'both')
                        tf_sync(condi,chanx,chany,:,:,1) = tmpsync;
                        tf_sync(condi,chanx,chany,:,:,2) = tmppli;
                    end
                    
                end % chany
            end % chanx
            if report_progress
                fprintf('done.\n')
            end
        end
    end % end condition loop
    
    %% db convert: condition-specific baseline
    if strcmp(baselinetype,'conspec')
        tf_pow = 10*log10( tf_pow ./ repmat(baselinedata,[ 1 1 1 length(times2save) ]) );
    elseif strcmp(baselinetype,'conavg')
        tf_pow = 10*log10( tf_pow ./ repmat(mean(baselinedata,1),[ nconds 1 1 length(times2save) ]) );
    end
    
    %% save results   
    if save_output
        chanlocs=ALLEEG(1).chanlocs;
        n=[ALLEEG.trials];
        if ~isempty(seeds)
            save(outputfilename,'tf_pow','tf_phase','tf_sync','frex','times2save','chanlocs','n','seeds');
        else
            save(outputfilename,'tf_pow','tf_phase','frex','times2save','chanlocs','n','seeds');            
        end
    end
    
    %% plot output
    if ~isempty(plot_output)
        
        tfwin = [plot_output.time; plot_output.freq];
        t2plot=dsearchn(times2save',tfwin(1,:)')';
        f2plot=dsearchn(frex',tfwin(2,:)')';
        chan2plot=[];
        for ch=1:length(plot_output.chan)
            chan2plot(ch) = find(strcmpi({EEG.chanlocs.labels},plot_output.chan(ch)));
        end

        figure('position',[10 10 400 600])
        subplot(421)
        contourf(times2save,frex,squeeze(mean(mean(tf_pow(:,chan2plot,:,:),1),2)),50,'linecolor','none')
        cl=get(gca,'clim');
        set(gca,'yscale',scale,'ytick',round(frex(1:4:end)),'clim',[-1*max(abs(cl)) max(abs(cl))])
        colorbar
        rectangle('position',[tfwin(1,1) tfwin(2,1) tfwin(1,2)-tfwin(1,1)  tfwin(2,2)-tfwin(2,1)])
        title([EEG.chanlocs(chan2plot).labels ' conavg pow'])
        
        subplot(4,2,[2 4])
        topoplot(squeeze(mean(mean(mean( tf_pow(:,:,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),1),3),4)),EEG.chanlocs,'electrodes','off','emarker2',{chan2plot(:),'o','k',5,1},'whitebk','on','maplimits',[-1*max(abs(cl)) max(abs(cl))]);
        title([num2str(tfwin(1,1)) '-' num2str(tfwin(1,2)) ' ms; ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz pow'])

        subplot(423)
        plot(times2save,squeeze(mean(mean( tf_pow(:,chan2plot,f2plot(1):f2plot(2),:),2),3)))
        set(gca,'xlim',[times2save(1) times2save(length(times2save))])
        title([EEG.chanlocs(chan2plot).labels ' ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz pow'])

        subplot(425)
        contourf(times2save,frex,squeeze(mean(mean(tf_phase(:,chan2plot,:,:),1),2)),50,'linecolor','none')
        cl=get(gca,'clim');
        set(gca,'yscale',scale,'ytick',round(frex(1:4:end)),'clim',[-1*max(abs(cl)) max(abs(cl))])
        colorbar
        rectangle('position',[tfwin(1,1) tfwin(2,1) tfwin(1,2)-tfwin(1,1)  tfwin(2,2)-tfwin(2,1)])
        title([EEG.chanlocs(chan2plot).labels ' conavg phase'])

        subplot(4,2,[6 8])
        title([num2str(tfwin(1,1)) '-' num2str(tfwin(1,2)) ' ms; ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz pow'])
        topoplot(squeeze(mean(mean(mean( tf_phase(:,:,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),1),3),4)),EEG.chanlocs,'electrodes','off','emarker2',{chan2plot(:),'o','k',5,1},'whitebk','on','maplimits',[-1*max(abs(cl)) max(abs(cl))]);
        subplot(427)
        plot(times2save,squeeze(mean(mean( tf_phase(:,chan2plot,f2plot(1):f2plot(2),:),2),3)))
        set(gca,'xlim',[times2save(1) times2save(length(times2save))])
        title([EEG.chanlocs(chan2plot).labels ' ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz phase'])
%%
        if strcmp(connectivity,'pli') || strcmp(connectivity,'ispc')
            
            tmpsync = tf_sync;
            tmpsync(isnan(tmpsync))=0;
                        
            % just plot from the first seed
            figure('position',[10 10 400 600])
            subplot(421)
            contourf(times2save,frex,squeeze(mean(mean(tmpsync(:,1,chan2plot,:,:),1),3)),50,'linecolor','none')
            set(gca,'yscale',scale,'ytick',round(frex(1:4:end)))
            colorbar
            rectangle('position',[tfwin(1,1) tfwin(2,1) tfwin(1,2)-tfwin(1,1)  tfwin(2,2)-tfwin(2,1)])
            title([[EEG.chanlocs(chan2plot).labels] '-' seeds{1}])
            
            subplot(4,2,[2 4])
            if strcmp(connectivity,'ispc')
                baset = dsearchn(times2save',[-200;-100])';
                actsync = squeeze(mean(mean(mean( tmpsync(:,1,:,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),1),4),5));
                basesync = squeeze(mean(mean(mean( tmpsync(:,1,:,f2plot(1):f2plot(2),baset(1):baset(2)),1),4),5));
                topoplot(actsync-basesync,EEG.chanlocs,'electrodes','off','emarker2',{[chan2plot(:)' find(strcmpi({EEG.chanlocs.labels},seeds(1)))],'o','k',5,1},'whitebk','on');
            else
                topoplot(squeeze(mean(mean(mean( tmpsync(:,1,:,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),1),4),5)),EEG.chanlocs,'electrodes','off','emarker2',{[chan2plot(:)' find(strcmpi({EEG.chanlocs.labels},seeds(1)))],'o','k',5,1},'whitebk','on');
            end
            
            title([num2str(tfwin(1,1)) '-' num2str(tfwin(1,2)) ' ms; ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz sync'])
            
            subplot(423)
            plot(times2save,squeeze(mean(mean( tmpsync(:,1,chan2plot,f2plot(1):f2plot(2),:),3),4)))
            set(gca,'xlim',[times2save(1) times2save(length(times2save))])
            title([num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz sync'])
            
        elseif strcmp(connectivity,'both')
            
            tmpsync = tf_sync;
            tmpsync(isnan(tmpsync))=0;
            
            % just plot from the first seed
            figure('position',[10 10 400 600])
            subplot(421)
            contourf(times2save,frex,squeeze(mean(mean(tmpsync(:,1,chan2plot,:,:,1),1),3)),50,'linecolor','none')
            set(gca,'yscale',scale,'ytick',round(frex(1:4:end)))
            colorbar
            rectangle('position',[tfwin(1,1) tfwin(2,1) tfwin(1,2)-tfwin(1,1)  tfwin(2,2)-tfwin(2,1)])
            title([[EEG.chanlocs(chan2plot).labels] '-' seeds{1}])
            
            subplot(4,2,[2 4])
            baset = dsearchn(times2save',[-200;-100])';
            actsync = squeeze(mean(mean(mean( tmpsync(:,1,:,f2plot(1):f2plot(2),t2plot(1):t2plot(2),1),1),4),5));
            basesync = squeeze(mean(mean(mean( tmpsync(:,1,:,f2plot(1):f2plot(2),baset(1):baset(2),1),1),4),5));
            topoplot(actsync-basesync,EEG.chanlocs,'electrodes','off','emarker2',{[chan2plot(:)' find(strcmpi({EEG.chanlocs.labels},seeds(1)))],'o','k',5,1},'whitebk','on');
            title([num2str(tfwin(1,1)) '-' num2str(tfwin(1,2)) ' ms; ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz ispc'])
            
            subplot(423)
            plot(times2save,squeeze(mean(mean( tmpsync(:,1,chan2plot,f2plot(1):f2plot(2),:,1),3),4)))
            set(gca,'xlim',[times2save(1) times2save(length(times2save))])
            title([num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz ispc'])
           
            subplot(425)
            contourf(times2save,frex,squeeze(mean(mean(tmpsync(:,1,chan2plot,:,:,2),1),3)),50,'linecolor','none')
            set(gca,'yscale',scale,'ytick',round(frex(1:4:end)))
            colorbar
            rectangle('position',[tfwin(1,1) tfwin(2,1) tfwin(1,2)-tfwin(1,1)  tfwin(2,2)-tfwin(2,1)])
            title([[EEG.chanlocs(chan2plot).labels] '-' seeds{1}])
            
            subplot(4,2,[6 8])
            topoplot(squeeze(mean(mean(mean( tmpsync(:,1,:,f2plot(1):f2plot(2),t2plot(1):t2plot(2),2),1),4),5)),EEG.chanlocs,'electrodes','off','emarker2',{[chan2plot(:)' find(strcmpi({EEG.chanlocs.labels},seeds(1)))],'o','k',5,1},'whitebk','on');
            title([num2str(tfwin(1,1)) '-' num2str(tfwin(1,2)) ' ms; ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz ispc'])
            
            subplot(427)
            plot(times2save,squeeze(mean(mean( tmpsync(:,1,chan2plot,f2plot(1):f2plot(2),:,2),3),4)))
            set(gca,'xlim',[times2save(1) times2save(length(times2save))])
            title([num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz dwpli'])
        end
    end
    
end

%%
end