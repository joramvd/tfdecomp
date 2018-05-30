function [tf_cluststat] = tfclustperm_chan(cfg,tfdat,varargin)

% unpack the cfg
v2struct(cfg);

nsubjects = size(tfdat,1);
if ~isfield(cfg,'nperm')
    nperm = 2000;
end

if ~isempty(varargin)
    tfdat = tfdat - varargin{1};
end

if isfield(cfg,'mask');
    
    fprintf('Computing spatial cluster based on input mask...\n');
    % make sure the begin and end time is the same between mask structure
    % and data
    if sum(ismember(dim.times,mask.time))<length(dim.times)
        
        tfdat = tfdat(:,:,:,ismember(dim.times,mask.time));
        mask.tmapthresh = mask.tmapthresh(:,ismember(mask.time,dim.times));
        
        dim.times = dim.times(ismember(dim.times,mask.time));        
        mask.time = mask.time(ismember(mask.time,dim.times));
    end
    
    clusts = bwconncomp(mask.tmapthresh);
    nmasks = clusts.NumObjects;
    fprintf('%i tf-clusters found in mask... Looping over tf-clusters:\n',nmasks);
    
    if ~strcmp(mask.cfg_prev.toi,'all');
        toi = dsearchn(dim.times',mask.cfg_prev.toi')';toi=toi(1):toi(2);
    else
        toi=1:length(dim.times);
    end
    if ~strcmp(mask.cfg_prev.foi,'all');
        foi = dsearchn(dim.freqs',mask.cfg_prev.foi')';foi=foi(1):foi(2);
    else
        foi=1:length(dim.freqs);
    end
    
else
    fprintf('Computing spatial cluster based on one predefined time-frequency window...\n');
    nmasks = 1;
end


%% loop over different clusters (or only one tf-window)

for maski = 1:nmasks
    
    if isfield(cfg,'mask')
        fprintf('Performing %i permutations of win/clust %i:\n',nperm,nmasks);
        if isfield(mask.cfg_prev,'mask') || strcmp(mask.cfg_prev.avgoverfreq,'yes')
            % if we use a statistical mask for averaging over time/freq, we
            % need to check whether that mask was a result of another
            % statistical mask, which resulted in averaging over one
            % frequency band, so our cluster will only be in time, not freq
            % same holds for a mask of a statistical test over time
            tfdat_reshaped = reshape(squeeze(mean(tfdat(:,:,foi,toi),3)),[nsubjects length(dim.chans(1:64)) length(toi)] );
        else
            tfdat_reshaped = reshape(tfdat(:,:,foi,toi),[nsubjects length(dim.chans(1:64)) length(toi)*length(foi)] );
        end
        X = squeeze(mean(tfdat_reshaped(:,:,clusts.PixelIdxList{maski}),3));
    else
        toi = dsearchn(dim.times',time)';toi=toi(1):toi(2);
        foi = dsearchn(dim.freqs',freq)';foi=foi(1):foi(2);
        X = squeeze(mean(mean(tfdat(:,:,foi,toi),3),4));
    end
    
    [~,pvalst,~,tmp] = ttest(X);
    tmap = squeeze(tmp.tstat);
    realmean = squeeze(mean(X));
    [p_fdr, p_masked] = fdr( pvalst, pval);
    
    % uncorrected pixel-level threshold
    threshmean = realmean;
    tmapthresh = tmap;
    if ~isfield(cfg,'tail')
        tmapthresh(abs(tmap)<tinv(1-pval/2,nsubjects-1))=0;
        threshmean(abs(tmap)<tinv(1-pval/2,nsubjects-1))=0;
    elseif strcmp(cfg.tail,'left')
        tmapthresh(tmap>-1.*tinv(1-pval,nsubjects-1))=0;
        threshmean(tmap>-1.*tinv(1-pval,nsubjects-1))=0;
    elseif strcmp(cfg.tail,'right')
        tmapthresh(tmap<tinv(1-pval,nsubjects-1))=0;
        threshmean(tmap<tinv(1-pval,nsubjects-1))=0;
    end
    
    %% real_part
    
    for pos_neg = [0 1];
        
        whichclust=0; % cluster identifier index
        clust_lab = cell(1,1);
        
        if pos_neg
            sign_chans = find(tmapthresh>0);
        else
            sign_chans = find(tmapthresh<0);
        end
        
        % first we loop over each significant channel and determine whether its
        % neighbours are also significant
        % these are stored in separate clusters
        for chani=sign_chans
            
            whichclust=whichclust+1;
            clust_lab{whichclust} = {dim.chans(chani).labels};
            
            % now loop through its neigbours to check the size of the cluster
            neighb_chans=zeros(1,length(neighbours(chani).neighblabel));
            for neighbi=1:length(neighbours(chani).neighblabel)
                neighb_chans(neighbi) = find(strcmpi({dim.chans.labels},neighbours(chani).neighblabel(neighbi)));
            end
            sign_neighb_chans = neighb_chans(ismember(neighb_chans,sign_chans));
            
            for neighbi=sign_neighb_chans
                clust_lab{whichclust} = [clust_lab{whichclust} {dim.chans(neighbi).labels}];
            end
        end
        
        %% it is likely that there is overlap
        
        keep_rollin = true;
        toggle=false;
        clusti=1;
        while keep_rollin
            
            for clustj = 1:length(clust_lab)
                
                if (clustj==clusti) && (length(clust_lab)>1)
                    continue;
                elseif (clustj==clusti) && (length(clust_lab)==1)
                    toggle=true;
                    continue
                end
                
                if ~isempty(intersect(clust_lab{clusti},clust_lab{clustj}))
                    new_clust_lab =  unique([clust_lab{[clusti clustj]}]);
                    clust_lab([clusti clustj])=[];
                    clust_lab = [clust_lab {new_clust_lab}];
                    toggle=false;
                    break
                else
                    toggle=true;
                    continue
                end
                
            end
            if toggle, clusti=clusti+1; end
            if clusti>length(clust_lab), keep_rollin=false; end
            
        end
        
        %% get the sum of t-values of the clusters
        
        elecs_idx = cell(1,length(clust_lab));
        sum_of_tvals = zeros(1,length(clust_lab));
        for clusti=1:length(clust_lab)
            elecs_idx{clusti} = zeros(1,length(clust_lab{clusti}));
            for chani=1:length(elecs_idx{clusti})
                elecs_idx{clusti}(chani) = find(strcmpi({dim.chans.labels}, clust_lab{clusti}{chani}));
            end
            sum_of_tvals(clusti) = sum(tmapthresh(elecs_idx{clusti}));
        end
        
        if pos_neg
            pos_clusters.tval = sum_of_tvals;
            pos_clusters.label = clust_lab;
        else
            neg_clusters.tval = sum_of_tvals;
            neg_clusters.label = clust_lab;
        end
        
    end
    
    %% permutation_part
    cluster_pval = pval;
    
    permuted_tsum_clusters = zeros(1,nperm);
    
    for permi=1:nperm
        
        if mod(permi,100)==0, fprintf('..%i\n',permi); end
        clabels = logical(sign(randn(nsubjects,1))+1);
        temp_permute = X; % X is the data structure and is assumed to be a difference score
        temp_permute(clabels,:,:)=temp_permute(clabels,:,:)*-1;
        
        [~,~,~,tmp] = ttest(squeeze(temp_permute));
        
        faketmap = squeeze(tmp.tstat);
        faketmap(abs(faketmap)<tinv(1-pval/2,nsubjects-1))=0;
        if ~isfield(cfg,'tail')
            faketmap(abs(faketmap)<tinv(1-pval/2,nsubjects-1))=0;
        elseif strcmp(cfg.tail,'left')
            faketmap(faketmap>-1.*tinv(1-pval,nsubjects-1))=0;
        elseif strcmp(cfg.tail,'right')
            faketmap(faketmap<tinv(1-pval,nsubjects-1))=0;
        end
        
        %%
        for pos_neg = [0 1];
            
            whichclust=0; % cluster identifier index
            clust_lab = cell(1,1);
            
            if pos_neg
                sign_chans = find(faketmap>0);
            else
                sign_chans = find(faketmap<0);
            end
            
            % first we loop over each significant channel and determine whether its
            % neighbours are also significant
            % these are stored in separate clusters
            for chani=sign_chans
                
                whichclust=whichclust+1;
                clust_lab{whichclust} = {dim.chans(chani).labels};
                
                % now loop through its neigbours to check the size of the cluster
                neighb_chans=zeros(1,length(neighbours(chani).neighblabel));
                for neighbi=1:length(neighbours(chani).neighblabel)
                    neighb_chans(neighbi) = find(strcmpi({dim.chans.labels},neighbours(chani).neighblabel(neighbi)));
                end
                sign_neighb_chans = neighb_chans(ismember(neighb_chans,sign_chans));
                
                for neighbi=sign_neighb_chans
                    clust_lab{whichclust} = [clust_lab{whichclust} {dim.chans(neighbi).labels}];
                end
            end
            
            %% there may be overlap in permuted clusters, too
            
            keep_rollin = true;
            toggle=false;
            clusti=1;
            while keep_rollin
                
                for clustj = 1:length(clust_lab)
                    
                    if (clustj==clusti) && (length(clust_lab)>1)
                        continue;
                    elseif (clustj==clusti) && (length(clust_lab)==1)
                        toggle=true;
                        continue
                    end
                    
                    if ~isempty(intersect(clust_lab{clusti},clust_lab{clustj}))
                        new_clust_lab =  unique([clust_lab{[clusti clustj]}]);
                        clust_lab([clusti clustj])=[];
                        clust_lab = [clust_lab {new_clust_lab}];
                        toggle=false;
                        break
                    else
                        toggle=true;
                        continue
                    end
                    
                end
                if toggle, clusti=clusti+1; end
                if clusti>length(clust_lab), keep_rollin=false; end
                
            end
            
            %%
            elecs_idx = cell(1,length(clust_lab));
            sum_of_tvals = zeros(1,length(clust_lab));
            for clusti=1:length(clust_lab)
                elecs_idx{clusti} = zeros(1,length(clust_lab{clusti}));
                for chani=1:length(elecs_idx{clusti})
                    elecs_idx{clusti}(chani) = find(strcmpi({dim.chans.labels}, clust_lab{clusti}{chani}));
                end
                sum_of_tvals(clusti) = sum(faketmap(elecs_idx{clusti}));
            end
            
            if pos_neg
                fake_pos_clusters.tval = sum_of_tvals;
            else
                fake_neg_clusters.tval = sum_of_tvals;
            end
            
        end
        
        %%
        all_tvals = [fake_pos_clusters.tval fake_neg_clusters.tval];
        permuted_tsum_clusters(permi) = max(abs(all_tvals));
        
    end
    fprintf('done!\n');
    clust_threshold = prctile(permuted_tsum_clusters,100-cluster_pval*100);
    
    %% cluster-size_threshold

    for clusti=1:length(pos_clusters.tval)
        pos_clusters.pvals(clusti) = length(find(permuted_tsum_clusters>abs(pos_clusters.tval(clusti))))/nperm;
    end
    for clusti=1:length(neg_clusters.tval)
        neg_clusters.pvals(clusti) = length(find(permuted_tsum_clusters>abs(neg_clusters.tval(clusti))))/nperm;
    end
    clusts2remove = find(abs(pos_clusters.tval)<clust_threshold);
    pos_clusters.tval(clusts2remove)=[];
    pos_clusters.label(clusts2remove)=[];
    pos_clusters.pvals(clusts2remove)=[];
    
    clusts2remove = find(abs(neg_clusters.tval)<clust_threshold);
    neg_clusters.tval(clusts2remove)=[];
    neg_clusters.label(clusts2remove)=[];
    neg_clusters.pvals(clusts2remove)=[];
    
    % channel indices of clusters, for plotting
    elecs2plot_lab = [neg_clusters.label{:},pos_clusters.label{:}]; % left electrodes; contralateral to right cues
    elecs2plot_idx = zeros(1,length(elecs2plot_lab));
    for chani=1:length(elecs2plot_lab)
        elecs2plot_idx(chani) = find(strcmpi({dim.chans.labels}, elecs2plot_lab{chani}));
    end
    
    %% now plot
    if strcmp(plot_output,'yes')
        labs2plot = [pos_clusters.label{:}, neg_clusters.label{:}];
        labs2plot_idx = zeros(1,length(labs2plot));
        for chani=1:length(labs2plot)
            labs2plot_idx(chani) = find(strcmpi({dim.chans.labels},labs2plot(chani)));
        end
        figure
        subplot(221)
        title('thresholded t-map')
        topoplot(tmapthresh,dim.chans,'electrodes','on','whitebk','on');
        subplot(222)
        title('t-map FDR')
        topoplot(tmap,dim.chans,'electrodes','off','emarker2',{find(p_masked),'o','w',5,1},'whitebk','on');
        subplot(223)
        title('t-map cluster-corrected')
        topoplot(tmap,dim.chans,'electrodes','off','emarker2',{labs2plot_idx,'o','w',5,1},'whitebk','on');
        
    end
    if isfield(cfg,'cmap')
        if ischar(cmap)
            colormap({cmap})
        else
            colormap(cmap)
        end
    end
    
    %% output
    
    tf_cluststat.pos_clusters(maski) = pos_clusters;
    tf_cluststat.neg_clusters(maski) = neg_clusters;
    tf_cluststat.all_clusters(maski).idx = elecs2plot_idx;
    tf_cluststat.topodat(maski).tmap = tmap;
    tf_cluststat.topodat(maski).realmap = realmean;
    
end
end
