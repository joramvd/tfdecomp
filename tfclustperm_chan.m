function [tf_cluststat] = tfclustperm_chan(cfg,X)

% unpack the cfg
v2struct(cfg);

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
        clust_lab{whichclust} = {channels(chani).labels};
        
        % now loop through its neigbours to check the size of the cluster
        neighb_chans=zeros(1,length(neighbours(chani).neighblabel));
        for neighbi=1:length(neighbours(chani).neighblabel)
            neighb_chans(neighbi) = find(strcmpi({channels.labels},neighbours(chani).neighblabel(neighbi)));
        end
        sign_neighb_chans = neighb_chans(ismember(neighb_chans,sign_chans));
        
        for neighbi=sign_neighb_chans
            clust_lab{whichclust} = [clust_lab{whichclust} {channels(neighbi).labels}];
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
            elecs_idx{clusti}(chani) = find(strcmpi({channels.labels}, clust_lab{clusti}{chani}));
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
fprintf('Performing %i permutations:\n',n_permutes);
cluster_pval = voxel_pval;

permuted_tsum_clusters = zeros(1,n_permutes);

for permi=1:n_permutes
    
    if mod(permi,100)==0, fprintf('..%i\n',permi); end
    clabels = logical(sign(randn(nsubjects,1))+1);
    temp_permute = X; % X is the data structure and is assumed to be a difference score
    temp_permute(clabels,:,:)=temp_permute(clabels,:,:)*-1;
    
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
    
%     figure
%     topoplot(faketmap,channels)

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
            clust_lab{whichclust} = {channels(chani).labels};
            
            % now loop through its neigbours to check the size of the cluster
            neighb_chans=zeros(1,length(neighbours(chani).neighblabel));
            for neighbi=1:length(neighbours(chani).neighblabel)
                neighb_chans(neighbi) = find(strcmpi({channels.labels},neighbours(chani).neighblabel(neighbi)));
            end
            sign_neighb_chans = neighb_chans(ismember(neighb_chans,sign_chans));
            
            for neighbi=sign_neighb_chans
                clust_lab{whichclust} = [clust_lab{whichclust} {channels(neighbi).labels}];
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
                elecs_idx{clusti}(chani) = find(strcmpi({channels.labels}, clust_lab{clusti}{chani}));
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

clusts2remove = find(abs(pos_clusters.tval)<clust_threshold);
pos_clusters.tval(clusts2remove)=[];
pos_clusters.label(clusts2remove)=[];

clusts2remove = find(abs(neg_clusters.tval)<clust_threshold);
neg_clusters.tval(clusts2remove)=[];
neg_clusters.label(clusts2remove)=[];

%% now plot
if plotfig
    labs2plot = [pos_clusters.label{:}, neg_clusters.label{:}];
    labs2plot_idx = zeros(1,length(labs2plot));
    for chani=1:length(labs2plot)
        labs2plot_idx(chani) = find(strcmpi({channels.labels},labs2plot(chani)));
    end
    figure
    topoplot(realmean,channels,'electrodes','off','emarker2',{labs2plot_idx,'o','w',5,1},'whitebk','on');
end

%% output

tf_cluststat.pos_clusters = pos_clusters;
tf_cluststat.neg_clusters = neg_clusters;

end
