%% template for level-1 analyses
% mikexcohen@gmail.com

%% user definitions

% IO definitions
readdir = '/path/to/post-processed/files/';
writdir = '/path/to/results/folder/';
outputbasename = 'tf';

% optional RT-locking (Assumes response codes immediately follow time-0-marker.
% Note also that trials with no response are REMOVED from the data.)
RTlocking = false; % true or false

% seeded phase synchronization (set to empty ("{}") for no analyses)
electrodes4seeded_synch = { 'fcz';'f5' };

% wavelet parameters
min_freq =  2;
max_freq = 40;
num_frex = 25;
wavelet_cycle_range = [ 3 12 ];
frequencyscaling = 'log'; % log or linear

% time parameters
times2save = -300:25:1300;
baselinetimerange = [-300 -100];

% files to use in analysis
filz = dir([ readdir '*final.mat' ]);

%% loop through files 

for subno=1:length(filz)
    
    if RTlocking, mn='_RTlocked'; else mn='_stimlocked'; end
    outputfilename = [ writdir filz(subno).name(1:end-4) mn '.mat' ];
    
    % the following line will skip over the subject if the output file already exists; 
    % comment this line to overwrite existing files
    if exist(outputfilename,'file'), continue; end
    
    load([ readdir filz(subno).name ])
    
    % match ALLEEG to EEG if only one dataset
    if ~exist('ALLEEG','var'), ALLEEG=EEG; end;
    if ~exist('EEG','var'), EEG=ALLEEG(1); end;
    
    %% additional data processing goes here...
    
    
    
    %% optional RT-locking
    
    if RTlocking
        
        % startpoint is a cell of vectors containing the real trial onset times (needed for baselining)
        startpoint    = cell(length(ALLEEG),1);
        baselinepoint = cell(length(ALLEEG),1);
        
        for condi=1:length(ALLEEG)
            
            startpoint{condi}    = zeros(length(ALLEEG(condi).epoch),1);
            baselinepoint{condi} = zeros(length(ALLEEG(condi).epoch),1);
            removetrial          = zeros(length(ALLEEG(condi).epoch),1);
            
            for ei=1:length(ALLEEG(condi).epoch)
                
                try
                    % find RT position
                    [j,rtloc  ] = min(abs(ALLEEG(condi).times - ALLEEG(condi).epoch(ei).eventlatency{find(cell2mat(ALLEEG(condi).epoch(ei).eventlatency)==0)+1}));
                    [j,zeroloc] = min(abs(ALLEEG(condi).times - 0 ));
                    
                    % find start point and replace with RT-aligned data
                    startpoint{condi}(ei)    = rtloc-zeroloc+1;
                    baselinepoint{condi}(ei) = zeroloc - (rtloc-zeroloc);
                    ALLEEG(condi).data(:,1:end-startpoint{condi}(ei)+1,ei) = ALLEEG(condi).data(:,startpoint{condi}(ei):end,ei);
                catch me, removetrial(ei)=1;
                end
            end
            
            % now remove trials with no response
            if sum(removetrial)>0
                ALLEEG(condi) = pop_select(ALLEEG(condi),'notrial',find(removetrial));
            end
        end % end condition loop
    end
    
    %% setup wavelet parameters
    
    if frequencyscaling(2)=='o'
        frex=logspace(log10(min_freq),log10(max_freq),num_frex);
    elseif frequencyscaling(2)=='i'
        frex=linspace(min_freq,max_freq,num_frex);
    else
        error('Unknown frequency scaling option!')
    end
        
    
    % gaussian width and time
    s=logspace(log10(wavelet_cycle_range(1)),log10(wavelet_cycle_range(2)),num_frex)./(2*pi.*frex);
    t=-EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
    
    % fft and convolution details
    Ltapr  = zeros(size(ALLEEG));
    Ldata  = zeros(size(ALLEEG));
    Lconv1 = zeros(size(ALLEEG));
    Lconv  = zeros(size(ALLEEG));
    
    for condi=1:length(ALLEEG)
        Ltapr(condi)  =  ALLEEG(condi).pnts;
        Ldata(condi)  =  ALLEEG(condi).pnts*ALLEEG(condi).trials;
        Lconv1(condi) =  Ldata(condi)+Ltapr(condi)-1;
        Lconv(condi)  =  pow2(nextpow2( Lconv1(condi) ));
    end
    
    wavelets = zeros(num_frex,length(t));
    for fi=1:num_frex
        wavelets(fi,:)=exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*s(fi)^2));
    end
    
    % setup time indexing
    times2saveidx = zeros(size(times2save));
    for ti=1:length(times2save)
        [~,times2saveidx(ti)]=min(abs(EEG.times-times2save(ti)));
    end    
    
    % baseline time indices
    [j,basetimeidx(1)] = min(abs(EEG.times-baselinetimerange(1)));
    [j,basetimeidx(2)] = min(abs(EEG.times-baselinetimerange(2)));
    
    % initialize output and hidden-layer TF matrices
    tf           = zeros(length(ALLEEG),EEG.nbchan,num_frex,length(times2saveidx),2);
    seededsynch  = zeros(length(ALLEEG),length(electrodes4seeded_synch),EEG.nbchan,num_frex,length(times2saveidx));
    regcoefs     = zeros(length(ALLEEG),ALLEEG(1).nbchan,num_frex,length(times2save));
    
    baselinedata = zeros(length(ALLEEG),ALLEEG(1).nbchan,num_frex);
    allphasevals = cell(size(ALLEEG));
    
    %% Now decompose
    
    chantoc=tic;
    
    % loop around channels
    for chani=1:ALLEEG(1).nbchan
        
        % finally, get FFT of condition-specific data
        EEGfft = cell(length(ALLEEG),1);
        for condi=1:length(ALLEEG)
            EEGfft{condi} = fft(reshape(ALLEEG(condi).data(chani,:,:),1,ALLEEG(condi).pnts*ALLEEG(condi).trials),Lconv(condi));
        end
        
        %% run standard wavelet analysis for low frequencies
        
        for fi=1:num_frex
            for condi=1:length(ALLEEG)
                
                % convolve and get analytic signal
                m = ifft(EEGfft{condi}.*fft(wavelets(fi,:),Lconv(condi)),Lconv(condi));
                m = m(1:Lconv1(condi));
                m = reshape(m(floor((Ltapr(condi)-1)/2):end-1-ceil((Ltapr(condi)-1)/2)),ALLEEG(condi).pnts,ALLEEG(condi).trials);
                
                % enter into TF matrix
                tf(condi,chani,fi,:,1) = mean(abs(m(times2saveidx,:)).^2,2)';
                tf(condi,chani,fi,:,2) = abs(mean(exp(1i*angle(m(times2saveidx,:))),2))';
                
                % also need single-trial data for baselining
                if RTlocking
                    for ti=1:size(m,2)
                        baselinedata(condi,chani,fi) = baselinedata(condi,chani,fi) + mean(abs(m(baselinepoint{condi}(ti)-min(baselinepoint{condi}(ti)-1,basetimeidx(1)):baselinepoint{condi}(ti)-basetimeidx(2),ti).^2));
                    end
                    baselinedata(condi,chani,fi) = baselinedata(condi,chani,fi)./ALLEEG(condi).trials; % average across trials
                else
                    baselinedata(condi,chani,fi) = mean(mean(abs(m(basetimeidx(1):basetimeidx(2),:)).^2,1),2);
                end
                
                %% optional phase population for phase synchronization (computed later)
                
                if ~isempty(electrodes4seeded_synch)
                    
                    % initialize
                    if chani==1 && fi==1
                        allphasevals{condi} = zeros(ALLEEG(1).nbchan,num_frex,length(times2save),ALLEEG(condi).trials);
                    end
                    
                    % populate
                    allphasevals{condi}(chani,fi,:,:) = angle(m(times2saveidx,:));
                end
                
                %% code for power regression
                
%                 for ti=1:length(times2save)
%                     [tempR,stats]=robustfit(regressors{condi},log(abs(m(times2saveidx(ti),:)').^2));
%                     regcoefs(condi,chani,fi,ti) = tempR(2)./stats.se(2);
%                 end % end time-loop
                
            end % end condition loop
        end % end frequency loop
        
        %% db convert
        
        if size(tf,1)==1
            tf(:,chani,:,:,1) = 10*log10( squeeze(tf(:,chani,:,:,1)) ./ repmat(squeeze(mean(baselinedata(:,chani,:),1)),[ size(tf,1) size(tf,4) ]) );
        else
            tf(:,chani,:,:,1) = 10*log10( squeeze(tf(:,chani,:,:,1)) ./ repmat(squeeze(mean(baselinedata(:,chani,:),1))',[ size(tf,1) 1 size(tf,4) ]) );
        end
        
        %% inter-site phase synchrony.
        
        for chanx=1:length(electrodes4seeded_synch)
            for chany=1:ALLEEG(1).nbchan
                for condi=1:length(ALLEEG)
                    seededsynch(condi,chanx,chany,:,:) = squeeze(abs(mean(exp(1i*( allphasevals{condi}(strcmpi(electrodes4seeded_synch{chanx},{ALLEEG(1).chanlocs.labels}),:,:,:)-allphasevals{condi}(chany,:,:,:) )),4)));
                end % condi
            end % chanj
        end % chani
        
    end % end channel loop

    %% ERPs

    erps = zeros(length(ALLEEG),ALLEEG(1).nbchan,ALLEEG(1).pnts);
    for condi=1:length(ALLEEG)
        erps(condi,:,:) = squeeze(mean(ALLEEG(condi).data,3));
    end
    erp_time = ALLEEG(1).times;

    
    %% save results
    
    chanlocs=ALLEEG(1).chanlocs;
    condition_labels={ALLEEG.setname};
    n=[ALLEEG.trials];

    save(outputfilename,'tf','frex','times2save','chanlocs','condition_labels','n','electrodes4seeded_synch','seededsynch','erps','erp_time');
end

%%
