
function outcl = feature_clust(spikes)

% outcl = feature_clust(spikes)
%
% Simple clustering based on feature activations. Currently this uses the
% spikes.active_feature and spikes.hosd_filt_peaks (unthresholded activation).
%
% The first round of clustering is based on the support of
% spikes.active_feature, such that two spikes are clustered together if
% they have supra-threshold activation in the same features.
% Clusters are then merged if the correlation of the mean value of unthresholded 
% activations (spikes.hosd_filt_peaks) between clusters exceeds 0.8.
% Finally, clusters are split using isosplit5 within the subset of
% features that are above threshold for at least 10% of the spikes within a
% given cluster.

  nthresh = ceil(size(spikes.active_feature,1)/2)+1; %Force merging of clusters smaller than this
% nthresh = 5; %Force merging of clusters smaller than this
% corrthresh = .6; 
mergethresh = .8; %Correlation  must be above this value to merge clusters
support_threshold = .1; %Consider a feature relevant if more than this proportion of spikes are active.
do_splitting = true; 
if do_splitting
    if isempty(which('isosplit5'))
        pth = fileparts(which(mfilename));
        addpath(fullfile(pth,'isosplit5','matlab'))
    end
end
%%% Preliminary clustering is based on the support of each feature.
[fclust,~,unqi] = unique(fliplr(spikes.active_feature'>0),'rows');
fclust = fliplr(fclust);

if isfield(spikes,'hosd_filt_pca')
    usefld = 'hosd_filt_pca';
else
    usefld = 'hosd_filt_peak';
end
%  usefld = 'active_feature';

%%% Now merge clusters based on the correlations among mean activations

mx = 1;
% for thri = 1:length(nthresh)
while any(mx>mergethresh)
    mfeat = zeros(size(fclust,1),size(spikes.active_feature,1));
    fsup= zeros(size(fclust,1),size(spikes.active_feature,1));
    for k = 1:size(fclust,1)
%               mfeat(k,:) = mean(spikes.active_feature(:,unqi==k),2);
            mfeat(k,:) = mean(spikes.(usefld)(:,unqi==k),2);
            fsup(k,:) = mean(spikes.active_feature(:,unqi==k)>0,2);
    end
     h = hist(unqi,1:max(unqi));
     nthresh = min(nthresh,max(h));
%     retain = h>=nthresh(thri);
    C = corr(mfeat');
    fretain = find(h>nthresh);
    trC = triu(C,1);
    [mx,mxi] = max(trC(fretain,:),[],1);
    mxi = fretain(mxi);
    retain = mx < mergethresh & h > nthresh;
%     [mx,mxi] = max(C(retain,:),[],1);
%     retain = retain | mx < corrthresh;
%     [mx,mxi] = max(C(retain,:),[],1);
    retain(mxi(~retain))=1;
    fretain = find(retain);
    [mx2,mxi2] = max(C(fretain,:),[],1);

    unqi = fretain(mxi2(unqi));
    [unq,~,unqi] = unique(unqi);
    fclust = fclust(unq,:);
    %params = spikes.params;
end
fclust = fsup > support_threshold;
%     retain = find(retain);
if do_splitting
    %%% Now clusters will be split using isosplit, based on a subset of
    %%% active features
    cl = zeros(size(unqi));
    subcl = zeros(size(unqi));
    nf = size(fclust,2);
    fprintf('\nClustering...')
    for k = 1:size(fclust,1)
         x = full(spikes.(usefld)(fclust(k,:),unqi==k))';
         if isempty(x)
             continue
         end
%         x = full(spikes.active_feature(fclust(k,:),unqi==k))';
    
        iso = isosplit5(x');
         
        cl(unqi==k) = 2^nf*iso + k;
        subcl(unqi==k) = iso;
        
    end
else
    cl = unqi;
end

[clusts,~,outcl] = unique(cl);


fprintf('%i clusters',length(clusts))
