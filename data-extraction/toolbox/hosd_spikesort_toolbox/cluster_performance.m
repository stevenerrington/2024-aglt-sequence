

function [out,target_cluster,spike_cluster] = cluster_performance(spike_cluster_in,target_clusters_in,Fs)

%Script to compute basic ground-truth comparisons

if nargin < 3 || isempty(Fs)
    Fs = spike_cluster_in.hos(1).sampling_rate;
end

if isnumeric(spike_cluster_in)
    spike_cluster_in = {spike_cluster_in};
end
if isnumeric(target_clusters_in)
    target_clusters_in = {target_clusters_in};
end

if iscell(spike_cluster_in)
    for k = 1:length(spike_cluster_in)
       if all(spike_cluster_in{k}-round(spike_cluster_in{k})==0)
           a = Fs;
       else
           a = 1;
       end
       spike_cluster.spike_sep(k).spike_times = double(spike_cluster_in{k})/a; %#ok<*AGROW>
       spike_cluster.spike_sep(k).cl = ones(size(spike_cluster_in{k}))*k;
       spike_cluster.spike_sep(k).cluster_count = length(spike_cluster_in{k});
    end
    spike_cluster.cl = [spike_cluster.spike_sep.cl];
    spike_cluster.spike_times = [spike_cluster.spike_sep.spike_times];    
    spike_cluster.Nclust = length(spike_cluster_in);
    spike_cluster.params = HOSD_default_params;
else
    spike_cluster = spike_cluster_in;
end

if iscell(target_clusters_in)
    for k = 1:length(target_clusters_in)
       if all(target_clusters_in{k}-round(target_clusters_in{k})==0)
           a = Fs;
       else
           a = 1;
       end
       target_clusters(k).spike_times = double(target_clusters_in{k})/a; %#ok<*AGROW>
    end
    
else
    target_clusters = target_clusters_in;
end
params = spike_cluster.params;
% [Ttru,tt] = chopper([-1 1]*spike_cluster.params.tolerance,[target_clusters.spike_times],Fs,length(ampfilt));%-mean(ttav),Fs);


maxn = round((max([spike_cluster.spike_times,target_clusters.spike_times])+.05)*Fs);


%  [srt,srti] = sort(spike_cluster.cl);
 
[Tdet,tt] = chopper([-1 1]*spike_cluster.params.tolerance,spike_cluster.spike_times,Fs,maxn);

view = length(unique(Tdet(:)))/maxn;
%  Tdet = mod(Tdet-1,maxn)+1;
qdet=zeros(maxn,1);
qdet(round(spike_cluster.spike_times*Fs)+1) = spike_cluster.cl;
detindx = zeros(maxn,1);
detindx(round(spike_cluster.spike_times*Fs)+1)=(1:length(spike_cluster.spike_times));
truindx = zeros(maxn,1);
truindx(round([target_clusters.spike_times]*Fs+1))=(1:length([target_clusters.spike_times]));

% qtru = zeros(maxn,1);
% qtru(round([target.spike_times]*Fs)) =1;

%  dbq = dbt(detrend([qdet',qtru'],0),Fs,10,'centerDC',false,'upsampleFx',4);
%  xc = ifft(mean(dbq.blrep(:,:,1).*conj(dbq.blrep(:,:,2))));
 
% [Ttru,tt] = chopper([-1 1]*params.tolerance,[target_clusters.spike_times],Fs,maxn);%-mean(ttav),Fs);
    qtru = zeros(maxn,1);

 Ttrus = {};
 for k = 1:length(target_clusters)
    [Ttru,tt] = chopper([-1 1]*spike_cluster.params.tolerance,[target_clusters(k).spike_times],Fs,maxn);%-mean(ttav),Fs);
    
    evtt = tt.*(qdet(Ttru)>0) + 0./qdet(Ttru);
    [~,mni] = min(abs(evtt));
    
     qtru(round(target_clusters(k).spike_times*Fs)+1) = k;

     evtt2 = tt.*(qtru(Ttru)>0) + 0./qtru(Ttru);
    [~,mni2] = min(abs(evtt2));
    
    target_clusters(k).clust = [0 1]'*sum(((1:length(tt))'==mni).*qdet(Ttru)); 
    target_clusters(k).clust(1,:) = k;
    target_clusters(k).detindx = sum(((1:length(tt))'==mni).*detindx(Ttru));
    target_clusters(k).truindx = sum(((1:length(tt))'==mni2).*truindx(Ttru));
%     target_clusters(k).truindx = sum((abs(tt)==min(abs(tt))).*truindx(Ttru));
  
     Ttrus{end+1}=Ttru;
 end
%  cl = [target_clusters.clust];
%  Q = (qtru(Tdet)>0).*hann(length(tt)).^4;
% ttav = (tt'*Q)./sum(Q);
% qt = sum(qtru(Tdet).*Q)./sum(Q);
Tdets = {};
 for k = 1:length(spike_cluster.spike_sep)
    Tdet = chopper([-1 1]*params.tolerance,[ spike_cluster.spike_sep(k).spike_times],Fs,maxn);
    
    evtt = tt.*(qtru(Tdet)>0);
    [mnt,mni] = min(abs(evtt) + 0./(qtru(Tdet)>0));
    evtt2 = tt.*(qdet(Tdet)>0);
    [~,mni2] = min(abs(evtt2) + 0./(qdet(Tdet)>0));
%     [~,mni2] = min(abs(tt));
   
    spike_cluster.spike_sep(k).clust = sum(((1:length(tt))'==mni).*qtru(Tdet)); 
    spike_cluster.spike_sep(k).clust(2,:) = k;
    spike_cluster.spike_sep(k).truindx = sum(((1:length(tt))'==mni).*truindx(Tdet));
    spike_cluster.spike_sep(k).detindx = sum(((1:length(tt))'==mni2).*detindx(Tdet));
%     spike_cluster.spike_sep(k).detindx = detindx(Tdet((1:length(tt))'==mni2,:))';
    Tdets{end+1} = Tdet;
 end
% cl = [ spike_cluster.spike_sep.clust];
cl = [target_clusters.clust,  spike_cluster.spike_sep.clust];
inds = [[target_clusters.detindx;target_clusters.truindx],[spike_cluster.spike_sep.detindx;spike_cluster.spike_sep.truindx]];

%%% Make sure there is no double counting.
%%%First, make sure not to discard any adc events.
[unq,unqi] = unique(inds([2 1],:)','rows');
cl = cl(:,unqi);
inds = inds(:,unqi);


[unq1,~,unqi1] = unique(inds(1,:)');
mult1 = find(crosstab(unqi1)>1 & unq1~=0);
multi1=find(ismember(unqi1,mult1));

[unq2,~,unqi2] = unique(inds(2,:)');
mult2 = find(crosstab(unqi2)>1 & unq2~=0);
multi=find(ismember(unqi1,mult1) | ismember(unqi2,mult2));
minds = inds(:,multi);

dscminds = find(diff([0,minds(2,:)])==0|minds(2,:)==0);
inds(:,multi(dscminds)) = [];
cl(:,multi(dscminds)) = [];

[a,b] = ndgrid(0:length(target_clusters),0:spike_cluster.Nclust);
clpad = [cl,[a(:)';b(:)']]; %Padding to make sure empty rows and coluns are not omitted in the confusion matrix
confMtx = crosstab(clpad(1,:),clpad(2,:))-1; %Effect of padding is removed by subtracting 1

 out.confMtx = confMtx;
 out.targetN = sum(sum(confMtx(2:end,1:end)));
 out.detectN = [1 ones(1,length(target_clusters))]*confMtx(:,2:end);
 out.all.Nmatch = sum(sum(confMtx(2:end,2:end)));
 out.all.Ntrue =sum( sum(confMtx(2:end,:)));
 out.all.Ndetect = sum(sum(confMtx(1:end,2:end)));
 out.all.NfalsePos = out.all.Ndetect - out.all.Nmatch;
 out.all.Nmiss = out.all.Ntrue - out.all.Nmatch;
 out.all.Accuracy = out.all.Nmatch./(out.all.Nmatch + out.all.NfalsePos + out.all.Nmiss); % Hit rate ignoring clustering
 out.all.Recall =out.all.Nmatch./out.all.Ntrue; % Hit rate ignoring clustering
 out.all.Precision = out.all.Nmatch./out.all.Ndetect; %Precision ignoring clustering
 out.all.note = 'Performance considering only detections (i.e. pooling all true and sorted units into 1 cluster, respectively)'; 
 %%% hit rate and precision counting the only the cluster with the highest
 %%% association as correct detection and everything else as
 %%% misclassification
%  [mx,mxi] = max(confMtx(2:end,2:end),[],2);
%  out.clusterhitrate = sum( confMtx(2:end,2:end)./sum(confMtx(2:end,1:end),2) .*((1:size(confMtx,2)-1)==mxi),2); %Hit rate within the most strongly associated cluster
%  out.cluster.Precision  = sum( confMtx(2:end,2:end)./sum(confMtx(1:end,2:end),1) .*((1:size(confMtx,2)-1)==mxi),2); % Precision for the same cluster
 % Computing accuracy following Magland 2020 (eqs 1-3)
Ntrue = sum(confMtx(2:end,1:end),2);
Nmatch = confMtx(2:end,2:end);
Ndetect = sum(confMtx(1:end,2:end));
Nmiss = Ntrue-Nmatch;
Nfalsepos = Ndetect - Nmatch;
Accuracy = Nmatch./(Nmatch + Nmiss + Nfalsepos);
Precision = Nmatch./Ndetect;
Recall = Nmatch./Ntrue;
force_one2one = false;
[row,col] = match_clusters(Accuracy,[],[]); %#ok<*UNRCH>
ncl = min(length(row),length(col));
one2one = zeros(size(row'));
one2one(row(1:ncl)) = col(1:ncl);
[~,mxAi] = max(Accuracy,[],2);
if force_one2one
    %%% Rather than using the max value, this enforces a 1-to-1 mapping between
    %%% sorted and true clusters by aligning accuracy  with the main diagonal
    %%% as best as possible.
    mxi = one2one;
else    
    mxi = mxAi;
end
%%% row and column vectors to sort Accuracy in descending diagonal order
out.AccuracySortingIndices = {row,col}; %For the sorted array, use AccSorted = out.full.Accuracy(out.AccuracySortingIndices{:});
out.accSortSorted = col;
out.one2one_assignment = one2one;
out.max_assignment = mxAi;
out.used_assignment = mxi;
%%% Stats for every pairing of true- to sorted- to cluster
out.full.Accuracy = Accuracy;
out.full.Precision = Precision;
out.full.Recall = Recall;
out.full.Nmatch = Nmatch;
out.full.Nmiss = Nmiss;
out.full.NfalsePos = Nfalsepos;
out.full.note = 'Stats for every pairing of true- and sorted-cluster';

%%% Stats pairing each true cluster to the sorted cluster with the highest
%%% accuracy value
out.cluster.Nmatch = sum((b(2:end,2:end)==mxi).*Nmatch,2);
out.cluster.Nmiss = sum((b(2:end,2:end)==mxi).*Nmiss,2);
out.cluster.NfalsePos = sum((b(2:end,2:end)==mxi).*Nfalsepos,2);
out.cluster.Accuracy = sum((b(2:end,2:end)==mxi).*Accuracy,2);
out.cluster.Precision = sum((b(2:end,2:end)==mxi).*Precision,2);
out.cluster.Recall = sum((b(2:end,2:end)==mxi).*Recall,2);
out.cluster.note ='Stats pairing each true cluster to a sorted cluster based on accuracy. See .used_assignment for the assignment.'; 

p = Ntrue./sum(Ntrue);
%%%Simple average of cluster stats
out.av.Accuracy = mean(out.cluster.Accuracy);
out.av.Precision = mean(out.cluster.Precision);
out.av.Recall = mean(out.cluster.Recall);
out.av.note = 'The simple average of performance stats in .cluster';

%Average weighted by true cluster size
out.wgtav.Accuracy = out.cluster.Accuracy'*p;
out.wgtav.Precision = out.cluster.Precision'*p;
out.wgtav.Recall = out.cluster.Recall'*p;
out.wgtav.note = 'Weighted average of performance stats in .cluster, weighted by cluster size';

%Stats obtained by pooling the match, miss and false positive numbers
%across clusters.
out.pool.Accuracy = sum(out.cluster.Nmatch)./sum(out.cluster.Nmatch+out.cluster.Nmiss+out.cluster.NfalsePos);
out.pool.Precision = sum(out.cluster.Nmatch)./sum(out.cluster.Nmatch+out.cluster.NfalsePos);
out.pool.Recall = sum(out.cluster.Nmatch)./sum(out.cluster.Nmatch+out.cluster.Nmiss);
out.pool.note = 'Performance stats by pooling values for hits, misses and false positives across clusters.';
% 

 out.view = view;
 out.target_id=(0:length(target_clusters))';
 out.cluster_id = (0:spike_cluster.Nclust);
 
%  out.SNRamp = SNRamp;
%  out.SNRsd = SNRsd;
%  out.F1 = 2*out.hitrate.*out.precision./(out.hitrate+out.precision);
%  out.FMindx = sqrt(out.hitrate*out.precision);
%  out.clusterF1 = 2*out.clusterhitrate.*out.cluster.Precision./(out.clusterhitrate+out.cluster.Precision);
%  out.clusterFM = sqrt(out.clusterhitrate.*out.cluster.Precision);
%  out.classF1 = 2*out.classhitrate.*out.classprecision./(out.classhitrate+out.classprecision);
%  out.classFM = sqrt(out.classhitrate.*out.classprecision);
%  
if nargout >1
    target_cluster.spike_sep = target_clusters;
    target_cluster.q= qtru;
    target_cluster.T = [Ttrus{:}];
    
    spike_cluster.q= qdet;
    spike_cluster.T = [Tdets{:}];
end
