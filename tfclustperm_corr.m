function [tf_cluststat] = tfclustperm_corr(cfg,X,Y)
v2struct(cfg);
% a few defaults
if ~isfield(cfg,'pval')
    pval = 0.05;
end
if ~isfield(cfg,'nperm')
    nperm = 2000;
end

pixel_pval = pval;
cluster_pval = pval;
nfreqs = length(dim.freqs);
ntimes = length(dim.times);

if strcmp(foi,'all')
    X = reshape(X,[nsubjects nfreqs*ntimes]);
    Y = repmat(Y,[1 nfreqs*ntimes]);
end
if strcmp(avgoverfreq,'yes')
    fx=dsearchn(dim.freqs',foi')';fx=fx(1):fx(2);
    X = squeeze(mean(X(:,fx,:),2));
    Y = repmat(Y,[1 ntimes]);
end 

rho_x = (sum(Y.*X) - nsubjects*mean(Y).*mean(X)) ./ (sqrt(sum(Y.^2)-nsubjects*mean(Y).^2) .* sqrt(sum(X.^2)-nsubjects*mean(X).^2));
rho_t = rho_x.*sqrt((nsubjects-2)./(1-rho_x.^2)); % +/- Inf where rho == 1
rho_p = 2*tcdf(-abs(rho_t),nsubjects-2);

% threshold real data
rho_thresh = rho_t; rho_thresh(rho_p>pixel_pval)=0;
if strcmp(foi,'all')
    rho_thresh = reshape(rho_thresh,[nfreqs ntimes]);
end

real_clusts = bwconncomp(rho_thresh);
temp_r = rho_thresh(:);
real_tsum = zeros(1,real_clusts.NumObjects);
for clusti=1:real_clusts.NumObjects
    real_tsum(clusti) = sum(abs(temp_r(real_clusts.PixelIdxList{clusti}).*sqrt((nsubjects-2)./(1-temp_r(real_clusts.PixelIdxList{clusti}).^2))));
end

%% now permute
fprintf('Performing %i permutations:\n',nperm);

null_tsum = zeros(1,nperm);
for permi=1:nperm
    if mod(permi,100)==0, fprintf('..%i\n',permi); end

    clabels = logical(sign(randn(nsubjects,1))+1);
    fake_Y = Y;
    fake_Y(clabels,:) = fake_Y(clabels,:).*-1;
    
    fake_rho_x = (sum(fake_Y.*X) - nsubjects*mean(fake_Y).*mean(X)) ./ (sqrt(sum(fake_Y.^2)-nsubjects*mean(fake_Y).^2) .* sqrt(sum(X.^2)-nsubjects*mean(X).^2));
    fake_rho_t = fake_rho_x.*sqrt((nsubjects-2)./(1-fake_rho_x.^2)); % +/- Inf where rho == 1
    fake_rho_p = 2*tcdf(-abs(fake_rho_t),nsubjects-2);
    
    % threshold real data
    fake_rho_thresh = fake_rho_t; fake_rho_thresh(fake_rho_p>pixel_pval)=0;
    if strcmp(foi,'all')
        fake_rho_thresh = reshape(fake_rho_thresh,[nfreqs ntimes]);
    end
    
    null_clusts = bwconncomp(fake_rho_thresh);
    temp_r = fake_rho_thresh(:);
    temp_tsum = zeros(1,null_clusts.NumObjects);
    for clusti=1:null_clusts.NumObjects
        temp_tsum(clusti) = sum(abs(temp_r(null_clusts.PixelIdxList{clusti}).*sqrt((nsubjects-2)./(1-temp_r(null_clusts.PixelIdxList{clusti}).^2))));
    end
    
    if null_clusts.NumObjects>0, null_tsum(permi) = max(temp_tsum); end
end
fprintf('..Done!\n');

clust_threshold = prctile(null_tsum,100-cluster_pval*100);
clust_pvals = zeros(1,length(real_tsum));
for cp=1:length(real_tsum)
    clust_pvals(cp) = length(find(null_tsum>real_tsum(cp)))/nperm;
end

%% identify clusters to remove
whichclusters2remove = find(real_tsum<clust_threshold);
% cluster-size thresholding
for clusti=1:length(whichclusters2remove)
    rho_thresh(real_clusts.PixelIdxList{whichclusters2remove(clusti)})=0;
end

if strcmp(foi,'all')
    tf_cluststat.realmap = reshape(rho_x,[nfreqs ntimes]);
    tf_cluststat.tmap = reshape(rho_t,[nfreqs ntimes]);
else
    tf_cluststat.realmap = rho_x;
    tf_cluststat.tmap = rho_t;
end
tf_cluststat.tmapthresh = rho_thresh;
[tf_cluststat.pvals,idx] = sort(clust_pvals,2,'ascend'); % gives the cluster with lowest p-value first
tf_cluststat.clustinfo = real_tsum(idx);
tf_cluststat.cfg_prev = cfg;

end