function out = roughclust(spikes)

% clust = roughclust(spikes)
%
% Clustering based on components recovered by HOSD.
%
% Input:
%
%    spikes - output of HOSD_spike_detection.
%
% Output
%
%   clust - structure with the following main fields
%       .Nclust : number of clusters
%       .spike_times : spike times in sec.
%       .cl : cluster label for each detected spike
%       .params : parameter structure
%       .wint : [225×1 double]
%       .spike_sep : 1 x Nclust structure with data for each cluster in the 
%                   following fields:
%             .id : cluster label
%             .cluster_count : number of spikes in cluster
%             .indices : indices for spikes in cluster
%             .spike_times : spike times for spikes in cluster
%             .ach : autocorrelation histogram with 1 ms bins.
%             .active_feature : thresholded feature activations for each
%                               .spike
%             .mwave : mean waveform
%             .sdwave : waveform std. dev.
%             .spike_waves : all spike waves within cluster. 
%       .acht : autocorrelation histogram bin times.
%       .hos : mvhos object containing the detection filters and feature
%             waveforms
%       .normalization : normalization applied to each channel
%       .D : Pairwise Mahalanobis distances between clusters 
%       .d2centers : Mahalanobis distances to cluster mean for each spike
%       .x : feature activations
%       .Mx : mean feature activation for each cluster
%       .SS : covariance of feature activations for each cluster
%       .project : projection of feature activations to plot separated clusters
%
% See also HOSD_SPIKE_DETECTION and MVHOSD

%C. Kovach 2022

params = spikes.params;

if isfield(spikes,'hosd_filt_pca')
    usefield = 'hosd_filt_pca';
else
    usefield = 'hosd_filt_peak';
end
if ~isfield(spikes,'hosd_filt_peak')
    spikes.hosd_filt_peak = squeeze(max(spikes.all_features));
end

cl = feature_clust(spikes);
ncl = max(cl);
binsz = .001;
ups = 3;%Upsampling to avoid binning adjacent spikes as 1
a = 1./(spikes.srate*binsz);
len = max(round(spikes.spike_indices*a*ups));
q = zeros(len,1);

qt = ifftshift((0:length(q)-1) - floor(length(q)/2))*binsz/ups;
qtd = qt(1:ups:end);

sm = zeros(size(q));
sm(1:ups)=1;
sm = circshift(sm,-floor(ups/2));
sm = fft(sm);
ach_limit = .02;
fprintf('\nGetting stats for cluster: ')
nfp = 0;
for k = 1:max(cl)
    nfp = fprintf([repmat('\b',1,nfp),'%i of %i'],k,ncl)-nfp;
    q(:) = 0;
    spki = find(cl==k);
    spt = spikes.spike_indices(spki);
    cli = round(spt*a*ups);
    q(cli) = 1;
    ach = ifft(abs(fft(q)).^2);
    ach(1) = 0;
    ach = ifft(sm.*fft(ach));
    ach = ach(1:ups:end);
    ach = fftshift(ach(abs(qt(1:ups:end))<ach_limit));
    clust(k).id = k; %#ok<*AGROW>
    clust(k).cluster_count = length(spki);
    clust(k).indices = spki;

    clust(k).spike_times = spt/spikes.srate;
    clust(k).ach = ach;
    clust(k).active_feature = spikes.active_feature(:,cl==k);
    if ~isempty(spikes.spike_waves)
        clust(k).mwave = squeeze(mean(spikes.spike_waves(:,cl==k,:),2));
        clust(k).sdwave = squeeze(std(spikes.spike_waves(:,cl==k,:),[],2));
        clust(k).spike_waves = spikes.spike_waves(:,cl==k,:);
%     else
%         Tcl = T(:,cl==k);
% 
%         for kk = 1:size(mua,2)
%             X = mua(Tcl+(kk-1)*size(mua,1));
%             clust(k).mwave(:,kk) = mean(X,2);
%             clust(k).sdwave(:,kk) = std(X,[],2);
%         end
    end
    
end
  
out.Nclust = length(clust);

for k = 1:length(clust)
    out.cl(clust(k).indices)=k;
    out.Mx(k,:) = mean(spikes.(usefield)(:,clust(k).indices),2);
    out.SS(:,:,k) = cov(spikes.(usefield)(:,clust(k).indices)');
end
out.params = params;
out.x = spikes.(usefield)';
out.spike_times = spikes.spike_indices./spikes.srate;
out.used_components = ones(1,size(spikes.active_feature,1));
out.wint = spikes.tt;
out.spike_sep = clust;
out.hos = spikes.hos;
out.normalization = spikes.normalization;
[out.D,out.d2centers,~,out.project,out.discarded,out.SS,out.Mx] = clusterdists(out.x,out.cl,Inf);
out.acht = fftshift(qtd(abs(qtd)<ach_limit));
