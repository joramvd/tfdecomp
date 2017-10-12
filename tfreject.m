function [rejtrials] = tfreject(cfg,varargin)

%% unpack cfg and setup settings

v2struct(cfg)

path_list_cell = regexp(path,pathsep,'Split');
if (exist(eeglab_path,'dir') == 7) && ~any(ismember(eeglab_path,path_list_cell))
    % addpath(genpath(eeglab_path),'-begin');
    addpath(genpath(eeglab_path));
end

% if not specified, set fields to false
fields2check = {'plot_output','report_progress','overwrite'};
for fieldi = 1:length(fields2check)
    if isfield(cfg,fields2check(fieldi))
        eval([fields2check{fieldi} ' = ' fields2check{fieldi} ';']);
        if strcmp(eval(fields2check{fieldi}),'none')
            eval([fields2check{fieldi} ' = false;']);
        end
    else
        eval([fields2check{fieldi} ' = false;']);
    end
end

% files to use in analysis
if isempty(varargin{1})
    filz = dir([ readdir filesep filename ]);
    if isempty(filz)
        error('No files to be loaded...!')
    end
else
    EEG = varargin{1}; % this implies the second argument of tfdecomp, if specified, is an ALLEEG structure
    filz=1;
end

% frequencies
if strcmp(scale,'log')
    frex=logspace(log10(frequencies(1)),log10(frequencies(2)),frequencies(3));
elseif strcmp(scale,'lin')
    frex=linspace(frequencies(1),frequencies(2),frequencies(3));
end
num_freqs = frequencies(3);

% gaussian width and time
npoints = srate*sum(abs(epochtime));
s=logspace(log10(cycles(1)),log10(cycles(2)),num_freqs)./(2*pi.*frex);
t=-npoints/srate/2:1/srate:npoints/srate/2-1/srate;

wavelets = zeros(num_freqs,length(t));
for fi=1:num_freqs
    wavelets(fi,:)=exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*s(fi)^2));
end

% setup time indexing
times2saveidx = zeros(size(times2save));
for ti=1:length(times2save)
    [~,times2saveidx(ti)]=min(abs(EEG.times-times2save(ti)));
end

rawconv = zeros(length(channels),num_freqs,length(times2save),EEG.trials);
Lconv = pow2(nextpow2( npoints*EEG.trials + npoints-1 ));

%% decompose
reverseStr='';
% loop around channels
for chani=channels
    
    if report_progress
        % display progress
        msg = sprintf('Decomposing channel %i/%i...',  chani,length(channels));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    try
        EEGfft = fft(reshape(EEG.data(chani,:,:),1,npoints*EEG.trials),Lconv);
    catch me
        warning('More channels specified than present in data, rest is set to zero')
        continue
    end
    
    % loop around frequencies
    for fi=1:num_freqs
        
        % convolve and get analytic signal
        m = ifft(EEGfft.*fft(wavelets(fi,:),Lconv),Lconv);
        m = m(1:(npoints*EEG.trials + npoints-1));
        m = reshape(m(floor((npoints-1)/2):end-1-ceil((npoints-1)/2)),npoints,EEG.trials);
        
        % populate
        rawconv(chani,fi,:,:) = m(times2saveidx,:);
        
    end % end frequency loop
end % end channel loop

%%
% additional trial rejection procedure for power outlier trials
rawpow = squeeze(mean( reshape((abs(rawconv).^2),length(channels)*num_freqs*length(times2save),EEG.trials), 1));
rejtrials=find(rawpow>median(rawpow)+3*std(rawpow));
fprintf('Removed an additional %i trials.\n', length(rejtrials));


end