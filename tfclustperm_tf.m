function [tf_cluststat] = tfclustperm_tf(cfg,X)

% ingredients from cfg

% a few defaults
if ~isfield(cfg,'plot_output')
    cfg.plot_output = 'no';
end
if ~isfield(cfg,'pval')
    cfg.pval = 0.05;
end
if ~isfield(cfg,'conn')
    cfg.conn = 8;
end
if ~isfield(cfg,'ndim')
    cfg.ndim = 1;
end
if ~isfield(cfg,'nperm')
    cfg.nperm = 2000;
end
if ~isfield(cfg,'test_statistic')
    cfg.test_statistic = 'sum';
end
if ~isfield(cfg,'avgoverfreq')
    cfg.avgoverfreq = 'no';
end

tffrex = cfg.freq;
tftime = cfg.time;
nFreq = numel(tffrex);
nTime = numel(tftime);
nsubjects = cfg.nsubjects;
voxel_pval   = cfg.pval;
cluster_pval = cfg.pval;

if ndims(X)==4 && size(X,ndims(X))==2, X=squeeze(X(:,:,:,1)-X(:,:,:,2)); end

if ndims(X)==2, X=permute(X,[1 3 2]); end

% if requested, average over frequency domain to do time-domain test
if strcmp(cfg.avgoverfreq,'yes');
    X = mean(X,2);
    nFreq = 1;
end

% note: try to use 1000 or more permutations for real data
n_permutes = cfg.nperm;

%% cluster test based on z-values through permutation testing

if strcmp(cfg.type , 'Z')
    
    % initialize null hypothesis matrices
    permuted_vals    = zeros(n_permutes,nFreq,nTime);
    max_clust_info   = zeros(n_permutes,1);

    for permi=1:n_permutes

        clabels = logical(sign(randn(nsubjects,1))+1);
        tlabels = logical(sign(randn(nsubjects,1))+1);
        flabels = logical(sign(randn(nsubjects,1))+1);


        temp_permute = X; % X is the data structure and is assumed to be a difference score
        temp_permute(clabels,:,:)=temp_permute(clabels,:,:)*-1;
        if cfg.ndim>1,  
            cutLoc = 2 + randperm(nTime-4); cutLoc=cutLoc(1);
            temp_permute(tlabels,:,:)=temp_permute(tlabels,:,[cutLoc:end 1:cutLoc-1]); 
        end

        if cfg.ndim>2,
            cutLoc = 1 + randperm(nFreq-2); cutLoc=cutLoc(1);
            temp_permute(flabels,:,:)=temp_permute(flabels,[cutLoc:end 1:cutLoc-1],:); 
        end

        permuted_vals(permi,:,:) = mean(temp_permute);
    end

    realmean = squeeze(mean(X));
    zmap = squeeze((mean(X)-mean(permuted_vals)) ./ std(permuted_vals));

    %%
    threshmean = squeeze(mean(X));
    threshmean(abs(zmap)<norminv(1-voxel_pval/2))=0;

    %%

    fprintf('Performing %i permutations:\n',n_permutes);

    % this time, the cluster correction will be done on the permuted data, thus
    % making no assumptions about parameters for p-values
    for permi = 1:n_permutes

        if mod(permi,100)==0, fprintf('..%i\n',permi); end
        % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
        fakecorrsz = squeeze((permuted_vals(permi,:,:)-mean(permuted_vals,1)) ./ std(permuted_vals,[],1) );
        fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval/2))=0;

        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(fakecorrsz,cfg.conn);
        if strcmp(cfg.test_statistic,'count')
            max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
        elseif strcmp(cfg.test_statistic,'sum')
            tmp_clust_sum = zeros(1,clustinfo.NumObjects);
            for ii=1:clustinfo.NumObjects
                tmp_clust_sum(ii) = sum(abs(fakecorrsz(clustinfo.PixelIdxList{ii})));
            end
            if  clustinfo.NumObjects>0, max_clust_info(permi) = max(abs(tmp_clust_sum)); end
        else
            error('Absent or incorrect test statistic input!');
        end

    end
    fprintf('..Done!\n');


    % apply cluster-level corrected threshold
    zmapthresh = zmap;
    % uncorrected pixel-level threshold
    zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval/2))=0;
    % find islands and remove those smaller than cluster size threshold
    clustinfo = bwconncomp(zmapthresh,cfg.conn);
    if strcmp(cfg.test_statistic,'count')
        clust_info = cellfun(@numel,clustinfo.PixelIdxList); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
    elseif strcmp(cfg.test_statistic,'sum')
        clust_info = zeros(1,clustinfo.NumObjects);
        for ii=1:clustinfo.NumObjects
            clust_info(ii) = sum(abs(zmapthresh(clustinfo.PixelIdxList{ii})));
        end
    end
    clust_threshold = prctile(max_clust_info,100-cluster_pval*100);

    % identify clusters to remove
    whichclusters2remove = find(clust_info<clust_threshold);

    % compute p-n value for all clusters
    %clusts2retain = find(clust_info>clust_threshold);
    %if ~isempty(clusts2retain)
    clust_pvals = zeros(1,length(clust_info));
    clust_act = clust_pvals;
    for cp=1:length(clust_info)
        clust_pvals(cp) = length(find(max_clust_info>clust_info(cp)))/n_permutes;
        clust_act(cp) = sum(zmapthresh(clustinfo.PixelIdxList{cp}));
    end
    %end
    
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end


    %% Generate figure: time-freq

    if strcmp(cfg.plot_output,'yes') && strcmp(cfg.avgoverfreq,'no');
        figure
        subplot(221)
        contourf(tftime,tffrex,squeeze(mean(X)),40,'linecolor','none')
        cmax=max(abs(get(gca,'clim')));
        axis square
        set(gca,'clim',[-cmax cmax],'yscale','log','ytick',[round(logspace(log10(min(tffrex)),log10(max(tffrex)),5))] )
        title('power map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')

        subplot(222)
        contourf(tftime,tffrex,zmap,40,'linecolor','none')
        cmax=max(abs(get(gca,'clim')));
        axis square
        set(gca,'clim',[-cmax cmax],'yscale','log','ytick',[round(logspace(log10(min(tffrex)),log10(max(tffrex)),5))] )
        title('unthresholded Z map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')

        subplot(223)
        contourf(tftime,tffrex,threshmean,40,'linecolor','none')
        cmax=max(abs(get(gca,'clim')));
        axis square
        set(gca,'clim',[-cmax cmax],'yscale','log','ytick',[round(logspace(log10(min(tffrex)),log10(max(tffrex)),5))] )
        title('Uncorrected power map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')

        subplot(224)
        contourf(tftime,tffrex,zmapthresh,40,'linecolor','none')
        cmax=max(abs(get(gca,'clim')));
        axis square
        set(gca,'clim',[-cmax cmax],'yscale','log','ytick',[round(logspace(log10(min(tffrex)),log10(max(tffrex)),5))] )
        title('Cluster-corrected Z map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    end
    %% Generate figure: time

    if strcmp(cfg.plot_output,'yes') && strcmp(cfg.avgoverfreq,'yes');
        figure
        subplot(121)
        plot(tftime,squeeze(mean(X)))
        axis square
        title('power map')
        xlabel('Time (ms)'), ylabel('Raw activity')

        hold on
        plotclust = bwconncomp(threshmean,cfg.conn);
        for blob=1:plotclust.NumObjects;
            plot(tftime(plotclust.PixelIdxList{blob}),threshmean(plotclust.PixelIdxList{blob}),'k','linewidth',4);
        end

        subplot(122)
        plot(tftime,zmap)
        axis square
        title('Z map')
        xlabel('Time (ms)'), ylabel('Z-val')

        hold on
        plotclust = bwconncomp(zmapthresh,cfg.conn);
        for blob=1:plotclust.NumObjects;
            plot(tftime(plotclust.PixelIdxList{blob}),zmapthresh(plotclust.PixelIdxList{blob}),'k','linewidth',4);
        end

    end

    %% output

    tf_cluststat.realmap = realmean;
    tf_cluststat.zmap = zmap;
    tf_cluststat.threshmean = threshmean;
    tf_cluststat.zmapthresh = zmapthresh;
    [tf_cluststat.pvals,idx] = sort(clust_pvals,2,'ascend'); % gives the cluster with lowest p-value first
    tf_cluststat.clustinfo = clust_act(idx);

    
%% cluster test based on t-tests

elseif strcmp(cfg.type, 'T');
    
    % initialize null hypothesis matrices
    %permuted_vals    = zeros(n_permutes,nsubjects,nFreq,nTime);
    max_clust_info   = zeros(n_permutes,1);
    
    %% real t-values
    
    [~,~,~,tmp] = ttest(X);
    tmap = squeeze(tmp.tstat);
    
    realmean = squeeze(mean(X));
    
    % uncorrected pixel-level threshold
    threshmean = realmean;
    tmapthresh = tmap;
    if ~isfield(cfg,'tail')
        tmapthresh(abs(tmap)<tinv(1-voxel_pval/2,nsubjects-1))=0;
        threshmean(abs(tmap)<tinv(1-voxel_pval/2,nsubjects-1))=0;
    elseif strcmp(cfg.tail,'left')
        tmapthresh(tmap>-1.*tinv(1-voxel_pval,nsubjects-1))=0;        
        threshmean(tmap>-1.*tinv(1-voxel_pval,nsubjects-1))=0;
    elseif strcmp(cfg.tail,'right')
        tmapthresh(tmap<tinv(1-voxel_pval,nsubjects-1))=0;        
        threshmean(tmap<tinv(1-voxel_pval,nsubjects-1))=0;
    end
    
    %%
    fprintf('Performing %i permutations:\n',n_permutes);
    
    for permi=1:n_permutes
        
        if mod(permi,100)==0, fprintf('..%i\n',permi); end
        
        clabels = logical(sign(randn(nsubjects,1))+1);
        tlabels = logical(sign(randn(nsubjects,1))+1);
        flabels = logical(sign(randn(nsubjects,1))+1);
        
        temp_permute = X; % X is the data structure and is assumed to be a difference score
        temp_permute(clabels,:,:)=temp_permute(clabels,:,:)*-1;
        if cfg.ndim>1,
            cutLoc = 2 + randperm(nTime-4); cutLoc=cutLoc(1);
            temp_permute(tlabels,:,:)=temp_permute(tlabels,:,[cutLoc:end 1:cutLoc-1]);
        end
        
        if cfg.ndim>2,
            cutLoc = 1 + randperm(nFreq-2); cutLoc=cutLoc(1);
            temp_permute(flabels,:,:)=temp_permute(flabels,[cutLoc:end 1:cutLoc-1],:);
        end
                
        %% permuted t-maps
        [~,~,~,tmp] = ttest(squeeze(temp_permute));
        
        faketmap = squeeze(tmp.tstat);
        faketmap(abs(faketmap)<tinv(1-voxel_pval/2,nsubjects-1))=0;
        if ~isfield(cfg,'tail')
            faketmap(abs(faketmap)<tinv(1-voxel_pval/2,nsubjects-1))=0;
        elseif strcmp(cfg.tail,'left')
            faketmap(faketmap>-1.*tinv(1-voxel_pval,nsubjects-1))=0;
        elseif strcmp(cfg.tail,'right')
            faketmap(faketmap<tinv(1-voxel_pval,nsubjects-1))=0;
        end

        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(faketmap,cfg.conn);
        if strcmp(cfg.test_statistic,'count')
            max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
        elseif strcmp(cfg.test_statistic,'sum')
            tmp_clust_sum = zeros(1,clustinfo.NumObjects);
            for ii=1:clustinfo.NumObjects
                tmp_clust_sum(ii) = sum(abs(faketmap(clustinfo.PixelIdxList{ii})));
            end
            if  clustinfo.NumObjects>0, max_clust_info(permi) = max(tmp_clust_sum); end
        else
            error('Absent or incorrect test statistic input!');
        end
        
    end
    fprintf('..Done!\n');
    
    %% apply cluster-level corrected threshold
    
    % find islands and remove those smaller than cluster size threshold
    clustinfo = bwconncomp(tmapthresh,cfg.conn);
    if strcmp(cfg.test_statistic,'count')
        clust_info = cellfun(@numel,clustinfo.PixelIdxList); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
    elseif strcmp(cfg.test_statistic,'sum')
        clust_info = zeros(1,clustinfo.NumObjects);
        for ii=1:clustinfo.NumObjects
            clust_info(ii) = sum(abs(tmapthresh(clustinfo.PixelIdxList{ii})));
        end
    end
    clust_threshold = prctile(max_clust_info,100-cluster_pval*100);
    
    % identify clusters to remove
    whichclusters2remove = find(clust_info<clust_threshold);
    
    % compute p-n value for all clusters
    %clusts2retain = find(clust_info>clust_threshold);
    %if ~isempty(clusts2retain)
    clust_pvals = zeros(1,length(clust_info));
    clust_act = clust_pvals;
    for cp=1:length(clust_info)
        clust_pvals(cp) = length(find(max_clust_info>clust_info(cp)))/n_permutes;
        clust_act(cp) = sum(tmapthresh(clustinfo.PixelIdxList{cp}));
    end
    %end
    
    % remove clusters
    for i=1:length(whichclusters2remove)
        tmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end
    

    %% Generate figure: time-freq
    
    if strcmp(cfg.plot_output,'yes') && strcmp(cfg.avgoverfreq,'no');
        figure
        subplot(221)
        contourf(tftime,tffrex,realmean,40,'linecolor','none')
        cmax=max(abs(get(gca,'clim')));
        axis square
        set(gca,'clim',[-cmax cmax],'yscale','log','ytick',[round(logspace(log10(min(tffrex)),log10(max(tffrex)),5))] )
        title('power map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')
        
        subplot(222)
        contourf(tftime,tffrex,squeeze(tmap),40,'linecolor','none')
        cmax=max(abs(get(gca,'clim')));
        axis square
        set(gca,'clim',[-cmax cmax],'yscale','log','ytick',[round(logspace(log10(min(tffrex)),log10(max(tffrex)),5))] )
        title('unthresholded t-map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')
        
        subplot(223)
        contourf(tftime,tffrex,threshmean,40,'linecolor','none')
        cmax=max(abs(get(gca,'clim')));
        axis square
        set(gca,'clim',[-cmax cmax],'yscale','log','ytick',[round(logspace(log10(min(tffrex)),log10(max(tffrex)),5))] )
        title('Uncorrected power map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')
        
        subplot(224)
        contourf(tftime,tffrex,tmapthresh,40,'linecolor','none')
        cmax=max(abs(get(gca,'clim')));
        axis square
        set(gca,'clim',[-cmax cmax],'yscale','log','ytick',[round(logspace(log10(min(tffrex)),log10(max(tffrex)),5))] )
        title('Cluster-corrected t-map')
        xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    end
    %% Generate figure: time
    
    if strcmp(cfg.plot_output,'yes') && strcmp(cfg.avgoverfreq,'yes');
        figure
        subplot(121)
        plot(tftime,squeeze(mean(X)))
        axis square
        title('power map')
        xlabel('Time (ms)'), ylabel('Raw activity')
        
        hold on
        plotclust = bwconncomp(threshmean,cfg.conn);
        for blob=1:plotclust.NumObjects;
            plot(tftime(plotclust.PixelIdxList{blob}),threshmean(plotclust.PixelIdxList{blob}),'k','linewidth',4);
        end
        
        subplot(122)
        plot(tftime,tmap)
        axis square
        title('T-stats')
        xlabel('Time (ms)'), ylabel('T-val')
        
        hold on
        plotclust = bwconncomp(tmapthresh,cfg.conn);
        for blob=1:plotclust.NumObjects;
            plot(tftime(plotclust.PixelIdxList{blob}),tmapthresh(plotclust.PixelIdxList{blob}),'k','linewidth',4);
        end
        
    end
    if isfield(cfg,'cmap')
        if ischar(cfg.cmap)
            colormap({cfg.cmap})
        else
            colormap(cfg.cmap)
        end
    end
    
    %% output
    
    tf_cluststat.realmap = realmean;
    tf_cluststat.tmap = tmap;
    tf_cluststat.threshmean = threshmean;
    tf_cluststat.tmapthresh = tmapthresh;
    [tf_cluststat.pvals,idx] = sort(clust_pvals,2,'ascend'); % gives the cluster with lowest p-value first
    tf_cluststat.clustinfo = clust_act(idx);
    
end
end