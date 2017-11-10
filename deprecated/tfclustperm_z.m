if strcmp(stat , 'Z')
    
    % initialize null hypothesis matrices
    permuted_vals    = zeros(nperm,nFreq,nTime);
    max_clust_info   = zeros(nperm,1);

    for permi=1:nperm

        clabels = logical(sign(randn(nSubjects,1))+1);
        tlabels = logical(sign(randn(nSubjects,1))+1);
        flabels = logical(sign(randn(nSubjects,1))+1);


        temp_permute = X; % X is the data structure and is assumed to be a difference score
        temp_permute(clabels,:,:)=temp_permute(clabels,:,:)*-1;
        if ndim>1,  
            cutLoc = 2 + randperm(nTime-4); cutLoc=cutLoc(1);
            temp_permute(tlabels,:,:)=temp_permute(tlabels,:,[cutLoc:end 1:cutLoc-1]); 
        end

        if ndim>2,
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

    fprintf('Performing %i permutations:\n',nperm);

    % this time, the cluster correction will be done on the permuted data, thus
    % making no assumptions about parameters for p-values
    for permi = 1:nperm

        if mod(permi,100)==0, fprintf('..%i\n',permi); end
        % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
        fakecorrsz = squeeze((permuted_vals(permi,:,:)-mean(permuted_vals,1)) ./ std(permuted_vals,[],1) );
        fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval/2))=0;

        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(fakecorrsz,conn);
        if strcmp(test_statistic,'count')
            max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
        elseif strcmp(test_statistic,'sum')
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
    clustinfo = bwconncomp(zmapthresh,conn);
    if strcmp(test_statistic,'count')
        clust_info = cellfun(@numel,clustinfo.PixelIdxList); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
    elseif strcmp(test_statistic,'sum')
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
        clust_pvals(cp) = length(find(max_clust_info>clust_info(cp)))/nperm;
        clust_act(cp) = sum(zmapthresh(clustinfo.PixelIdxList{cp}));
    end
    %end
    
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end


    %% Generate figure: time-freq

    if strcmp(plot_output,'yes') && strcmp(avgoverfreq,'no');
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

    if strcmp(plot_output,'yes') && strcmp(avgoverfreq,'yes');
        figure
        
        if ~isempty(varargin)
            
            subplot(121)
            [l,p] = boundedline(dim.times,squeeze(mean(X1)),squeeze(std(X1))./sqrt(nSubjects),'b','alpha','transparency',.1);
            outlinebounds(l,p);
            h(1)=l;
            [l,p] = boundedline(dim.times,squeeze(mean(X2)),squeeze(std(X1))./sqrt(nSubjects),'r','alpha','transparency',.1);
            outlinebounds(l,p);
            h(2)=l;
            legend(h,'conA','conB');
            yl=get(gca,'ylim');

            axis square
            title('Conditions')
            xlabel('Time (ms)'), ylabel('Raw activity')

            hold on
            plotclust = bwconncomp(threshmean,conn);
            for blob=1:plotclust.NumObjects;
                plot([tftime(plotclust.PixelIdxList{blob}(1)) tftime(plotclust.PixelIdxList{blob}(end))],[min(yl)+sum(abs(yl))/20 min(yl)+sum(abs(yl))/20],[.5 .5 .5],'linewidth',4);
            end

            subplot(122)
            [l,p] = boundedline(dim.times,squeeze(mean(X)),squeeze(std(X))./sqrt(nSubjects),'k','alpha','transparency',.1);
            outlinebounds(l,p);
            axis square
            title('Difference')
            xlabel('Time (ms)'), ylabel('T-val')

            hold on
            plotclust = bwconncomp(tmapthresh,conn);
            for blob=1:plotclust.NumObjects;
                plot([tftime(plotclust.PixelIdxList{blob}(1)) tftime(plotclust.PixelIdxList{blob}(end))],[min(yl)+sum(abs(yl))/20 min(yl)+sum(abs(yl))/20],'k','linewidth',4);
            end
        else
            [l,p] = boundedline(dim.times,squeeze(mean(X)),squeeze(std(X))./sqrt(nSubjects),'k','alpha','transparency',.1);
            outlinebounds(l,p);
            axis square
            title('Single condition')
            xlabel('Time (ms)'), ylabel('T-val')

            hold on
            plotclust = bwconncomp(threshmean,conn);
            for blob=1:plotclust.NumObjects;
                plot([tftime(plotclust.PixelIdxList{blob}(1)) tftime(plotclust.PixelIdxList{blob}(end))],[min(yl)+sum(abs(yl))/10 min(yl)+sum(abs(yl))/20],[.5 .5 .5],'linewidth',4);
            end
            hold on
            plotclust = bwconncomp(tmapthresh,conn);
            for blob=1:plotclust.NumObjects;
                plot([tftime(plotclust.PixelIdxList{blob}(1)) tftime(plotclust.PixelIdxList{blob}(end))],[min(yl)+sum(abs(yl))/20 min(yl)+sum(abs(yl))/20],'k','linewidth',4);
            end
            
        end
        
    end
    
    if isfield(cfg,'cmap')
        if ischar(cmap)
            cmap({cmap})
        else
            cmap(cmap)
        end
    end

    %% output

    tf_cluststat.realmap = realmean;
    tf_cluststat.subjmap = X;
    tf_cluststat.zmap = zmap;
    tf_cluststat.threshmean = threshmean;
    tf_cluststat.zmapthresh = zmapthresh;
    tf_cluststat.time = tftime;
    tf_cluststat.freq = tffrex;
    [tf_cluststat.pvals,idx] = sort(clust_pvals,2,'ascend'); % gives the cluster with lowest p-value first
    tf_cluststat.clustinfo = clust_act(idx);
    tf_cluststat.cfg_prev = cfg;

