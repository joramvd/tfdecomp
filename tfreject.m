function [rejtrials] = tfreject(cfg,eegdat)

%% unpack cfg and setup settings

v2struct(cfg)

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

outputfilename = [writdir filename];
disp(outputfilename)

% frequencies
if strcmp(scale,'log')
    frex=logspace(log10(frequencies(1)),log10(frequencies(2)),frequencies(3));
elseif strcmp(scale,'lin')
    frex=linspace(frequencies(1),frequencies(2),frequencies(3));
end
nfreqs = frequencies(3);

% gaussian width and time
ntimepoints = size(eegdat{1},2);
s=logspace(log10(cycles(1)),log10(cycles(2)),nfreqs)./(2*pi.*frex);
t=-ntimepoints/srate/2:1/srate:ntimepoints/srate/2-1/srate;

wavelets = zeros(nfreqs,length(t));
for fi=1:nfreqs
    wavelets(fi,:)=exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*s(fi)^2));
end

ntrials = cellfun(@(x) size(x,3),eegdat);
nconds = length(eegdat);

% setup time indexing
times2saveidx = zeros(size(times2save));
for ti=1:length(times2save)
    [~,times2saveidx(ti)]=min(abs(eegtime-times2save(ti)));
end

%% Now decompose
reverseStr='';

% loop around conditions
for condi=1:nconds
    
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
            
        end % end frequency loop
    end % end channel loop
    if report_progress
        fprintf('done.\n')
        reverseStr='';
    end
    
    % additional trial rejection procedure for power outlier trials
    % note: because we just average over all frequencies, timepoints and
    % channels, low-frequency will dominate due to 1/f; but this is in fact
    % what we want: high-frequency noise was already removed during
    % preprocessing (muscle artifacts); this step is meant to reject some
    % trials in which one channel showed e.g. a weird 'blip'; this is
    % likely reflected in all frequency bands
    % I tried some baseline-correction or other normalization procedures,
    % but the simple raw power exactly identifies the trials that can have
    % a dramatic impact on further analyses.
    
    rawpow = squeeze(mean( reshape(abs(rawconv).^2,length(channels)*nfreqs*length(times2save),ntrials(condi)), 1));
    rejtrials=(rawpow>median(rawpow)+ndev*std(rawpow));

    fprintf('Removed an additional %i trials.\n', length(rejtrials));
    
    save(outputfilename,'rejtrials');

end % end condition loop

%%


end