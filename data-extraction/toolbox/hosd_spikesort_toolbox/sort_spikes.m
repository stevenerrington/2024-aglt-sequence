function out = sort_spikes(input,params)

if ischar(input)
    input = load(input);

end

% pca_up = 1.5; %Make pca dims a multiple of the number of HOSD components
if isstruct(input)
    if nargin < 2
        params = input.params;
    end
    switch params.sorting_params.sort_on
        case 'filtered_waves' %This differs from hosd_filt_segments in that the filter is applied to the spike waveforms without deflation.
            ff = [input.hos.filterfft];
            ff = ff./sqrt(sum(sum(abs(ff).^2),3));
            ff = permute(ff,[1 4 3 2]);
            data = ifft(fft(input.spike_waves).*ff);
            data = fftshift(data(abs(input.hos(1).sampt/input.hos(1).sampling_rate)<2*params.tolerance,:,:,:),1);
            
        otherwise
            data = input.(params.sorting_params.sort_on);
    end
    dataorig= input.spike_waves;
else
    data = input;
    dataorig = data;
end

if ~isfield(params.sorting_params,'pca_up')
    params.sorting_params.pca_up = 1.5; %Make pca dims a multiple of the number of HOSD components
end
if size(data,3)>1 || size(data,4)>1
   % Concatenate across channels for svd 
   data = reshape(permute(data,[1 3 4 2]),size(data,1)*size(data,3)*size(data,4),size(data,2));
end

% [u,~] = svds(double(data'),ceil(params.keepc*pca_up));
[u,~] = svds(double(data'),ceil(sum(input.used_components)*params.sorting_params.pca_up));

%Decide what kind of spike sorting to do. 

if isa(params.sorting_params.algorithm,'function_handle') || params.sorting_params.use_evalclusters && ismember(params.sorting_params.algorithm,{'kmeans','linkage','gmdistribution','gmdist'})    
    clusterfun = params.sorting_params.algorithm;
else
    clusterfun = str2func(params.sorting_params.algorithm);
           
end

switch params.sorting_params.eval_criterion
   
    case 'chi2_criterion'
        
        eva = chi2_criterion(u,params.sorting_params.max_cluster_num,clusterfun,params.sorting_params.use_mog_estimates ,params.sorting_params.mog_ppweight ,[],[],params.sorting_params.cluster_number_penalty);
        out.clusterweight = eva.eval.gmdist.posterior(u);
        
    case {'bic_criterion','bic'}
        
        eva = bic_criterion(u,params.sorting_params.max_cluster_num,clusterfun,[],[],params.sorting_params.cluster_number_penalty);
        out.clusterweight = eva.eval.gmdist.posterior(u);
    case 'none'
        
        eva.OptimalY = clusterfun(u,params.sorting_params.max_cluster_num);
        eva.OptimalK = length(unique(eva.OptimalY));
    otherwise
        eva = evalclusters(u(:,1:params.keepc),clusterfun,params.sorting_params.eval_criterion, params.sorting_params.max_cluster_num,'klist',2:params.sorting_params.max_cluster_num);

   
        if eva.OptimalK == 2 %%% evalclust does not include 1 cluster in the comparison
            D = clusterdists(u(:,1:params.keepc),cl);
            if sqrt(D)<3
                eva.OptimalY(:) = 1;
                eva.OptimalK = 1;
            end
        end

end



%%% Retain the following additional fields from the input data structure, if they
%%% are present
checkfields = {'chan','block','study','date','time','subject','notes','channel'};
fldn= fieldnames(input);
getfields = fldn(ismember(lower(fldn),checkfields));
for k = 1:length(getfields)
    out.(getfields{k}) = input.(getfields{k});
end

out.Nclust = eva.OptimalK;
out.cl = eva.OptimalY;
out.clusteval = eva;
out.x=u;
out.params = params;


%%% Sort by size of the mean waveform
for k = 1:eva.OptimalK
   m  = mean(dataorig(:,out.cl==k,:),2);
   [~,mxi] = max(max(m)-min(m),[],3);
   M(:,k) = m(:,:,mxi);
  
end

if all(isnan(out.cl))
    out.cl = [];
    out.avg_waves = [];
    out.D = [];
    out.d2centers = [];
    out.project = [];
    out.discarded = [];
    out.SS = [];
    out.Mx = [];
    out.clusterweight = [];
    out.used_components = input.used_components;
    fprintf('\n\nThere are no spikes to cluster!');
else
    [~,srti] = sort(max(M)-min(M),'descend');
    isrt(srti) = 1:size(M,2) ;
    out.cl = isrt(out.cl);
    out.clusteval.OptimalY =isrt(out.clusteval.OptimalY) ;
    out.avg_waves = M(:,srti);
    [out.D,out.d2centers,~,out.project,out.discarded,out.SS,out.Mx] = clusterdists(u,out.cl,Inf);
    if isfield(out,'clusterweight')
        out.clusterweight = out.clusterweight(:,srti);
    end
end
out.used_components = input.used_components;

out.spike_times = (input.spike_indices-1)./input.srate; %First sample is t=0.
out.wint = input.tt;
for k = 1:eva.OptimalK
        spkt = out.spike_times(out.cl==k);
        isi = diff(spkt);
        out.spike_sep(k).cluster = k;
        out.spike_sep(k).inter_spike_qtiles(:,k) = quantile(isi,params.sorting_params.isi_qtiles);
%         out.spike_waves = input.spike_waves;
        out.spike_sep(k).cluster_count = sum(out.cl==k);
        out.spike_sep(k).spike_waves = input.spike_waves(:,out.cl==k,:);
        out.spike_sep(k).spike_times = out.spike_times(out.cl==k);
        out.spike_sep(k).active_feature=input.active_feature(:,out.cl==k);
end
out.hos = input.hos;
if isfield(input,'normalization')
    out.normalization = input.normalization;
end
if ~isempty(params.output_file.basename) || ~isempty(params.output_file.clusters)
    [pth,fn] = fileparts(params.output_file.clusters);
    
    if isempty(pth)
        pth = params.output_file.basepath;
    end
    if isempty(fn)
        fn = [params.output_file.basename,'_cluster.mat'];
    end
    out.params.output_file.clusters = fullfile(pth,fn);
    save(out.params.output_file.clusters,'-struct','out')
end
    
    
