function tfmultiplot(cfg,tfdata,dim)

v2struct(cfg);

%% check dimensions of data

if isfield(cfg,'single_subject')
    if single_subject
        tfdata = permute(tfdata,[length(size(tfdata))+1 1:length(size(tfdata))]);
    end
end
nsubjects = size(tfdata,1);

if (strcmp(metric,'pli') || strcmp(metric,'ispc')) && ndims(tfdata)<5
    error('connectivity plot requested but insufficient number of dimensions in data!')
end

%% reduce seed dimension in case of connectivity

if strcmp(metric,'pli') || strcmp(metric,'ispc')
    seed2plot = strcmpi(dim.cfg_prev.seeds,seed);
    tfdata = squeeze(nanmean(tfdata(:,:,seed2plot,:,:,:),3));
    if isfield(cfg,'single_subject') % squeezing removed the first singleton dimension again
        if single_subject
            tfdata = permute(tfdata,[length(size(tfdata))+1 1:length(size(tfdata))]);
        end
    end
end

%% metric labels

if strcmp(metric,'pow')
    mlab = 'Power (dB)';
elseif strcmp(metric,'phase')
    mlab = 'ITPC';
elseif strcmp(metric,'pli')
    mlab = 'dwPLI';
elseif strcmp(metric,'ispc')
    mlab = 'ISPC';
else
    error('Unknown metric requested...!')
end

%%
if isfield(cfg,'concomp')
    tfdata = mean(tfdata(:,concomp(1,:),:,:,:),2) - mean(tfdata(:,concomp(2,:),:,:,:),2);
    figtitle = 'Condition diff.';
else
    figtitle = 'Condition avg.';
end


%%

tfwin = [time; freq];
t2plot=dsearchn(dim.times',tfwin(1,:)')';
f2plot=dsearchn(dim.freqs',tfwin(2,:)')';
chan2plot=[];
for ch=1:length(chan)
    chan2plot(ch) = find(strcmpi({dim.chans.labels},chan(ch)));
end

if isfield(cfg,'chandiff')
    chan2plot2=[];
    for ch=1:length(chandiff)
        chan2plot2(ch) = find(strcmpi({dim.chans.labels},chandiff(ch)));
    end
    dat2plot = squeeze(mean(tfdata(:,:,chan2plot,:,:),3) - mean(tfdata(:,:,chan2plot2,:,:),3));
else
    dat2plot = squeeze(mean(tfdata(:,:,chan2plot,:,:),3));
end

if ndims(dat2plot)==3
    dat2plot = permute(dat2plot,[4 1 2 3]);
end
if ndims(dat2plot)==2
    dat2plot = permute(dat2plot,[3 4 1 2]);
end

figure('position',[50 50 600 400])
subplot(221)
contourf(dim.times,dim.freqs,squeeze(mean(mean(dat2plot,1),2)),50,'linecolor','none')
cl=get(gca,'clim');
set(gca,'yscale',scale,'ytick',round(dim.freqs(1:4:end)),'clim',[-1*max(abs(cl)) max(abs(cl))])
colorbar
rectangle('position',[tfwin(1,1) tfwin(2,1) tfwin(1,2)-tfwin(1,1)  tfwin(2,2)-tfwin(2,1)])
if isfield(cfg,'markevents')
    hold on
    for ev=1:length(markevents)
        plot([markevents(ev) markevents(ev)],[dim.freqs(1) dim.freqs(end)],'--k')
    end
end
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title([dim.chans(chan2plot).labels ' ' figtitle])

subplot(222)
topodat = tfdata;
topodat(isnan(topodat))=0;
topoplot(squeeze(mean(mean(mean(mean( topodat(:,:,:,f2plot(1):f2plot(2),t2plot(1):t2plot(2)),1),2),4),5)),dim.chans,'electrodes','off','emarker2',{chan2plot(:),'o','k',5,1},'whitebk','on','plotrad',0.65);
if isfield(cfg,'chandiff')
    topoplot([],dim.chans(chan2plot2),'style','blank','electrodes','on','plotrad',0.65,'emarker',{'o','w',5,1});
end
title([num2str(tfwin(1,1)) '-' num2str(tfwin(1,2)) ' ms; ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz'])

subplot(223)
if isfield(cfg,'linecols')
    hold on
    for li=1:size(linecols,1)
        if iscell(linecols)
            [l,p] = boundedline(dim.times,squeeze(mean(mean(dat2plot(:,li,f2plot(1):f2plot(2),:),1),3)),std(squeeze(mean(dat2plot(:,li,f2plot(1):f2plot(2),:,:),3)))./sqrt(nsubjects),linecols{li},'alpha','transparency',.1); hold on
        else
            [l,p] = boundedline(dim.times,squeeze(mean(mean(dat2plot(:,li,f2plot(1):f2plot(2),:),1),3)),std(squeeze(mean(dat2plot(:,li,f2plot(1):f2plot(2),:,:),3)))./sqrt(nsubjects),'alpha','transparency',.1); hold on
            l.Color = linecols(li,:);
            p.FaceColor = linecols(li,:);
        end
        h(li)=l;
        outlinebounds(l,p);
    end
else
    h=plot(dim.times,squeeze(mean(mean( dat2plot(:,:,f2plot(1):f2plot(2),:),1),3)));
end
yl = get(gca,'ylim');
set(gca,'xlim',[dim.times(1) dim.times(length(dim.times))])
if isfield(cfg,'markevents')
    hold on
    for ev=1:length(markevents)
        plot([markevents(ev) markevents(ev)],[yl(1) yl(2)],'--k')
    end
end
ylabel(mlab)
xlabel('Time (ms)')
title([dim.chans(chan2plot).labels ' ' num2str(tfwin(2,1)) '-' num2str(tfwin(2,2)) ' Hz'])
myLegend = legend(h,connames);

set(myLegend,'Units', 'pixels')
myOldLegendPos=get(myLegend,'Position');
hold on
h=subplot(224);
set(h,'Units', 'pixels')
set(h, 'Visible', 'off')
myPosition=get(h,'Position');
set(myLegend,'Position',[myPosition(1) myPosition(2) myOldLegendPos(3) myOldLegendPos(4)])

if isfield(cfg,'cmap')
    colormap(cmap)
end

    
end