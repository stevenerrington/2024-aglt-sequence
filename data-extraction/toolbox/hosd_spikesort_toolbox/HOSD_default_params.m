function params = HOSD_default_params(paramsin) 

%%%%
%%%% Below are default parameter values for preprocessing and HOSD-based detection %%%%
%%%%

params.windur = .015; % HOS estimator window duration (sec)
params.lowpass = 4000; % lowpass cutoff ( Hz )
params.highpass = 200; % highpass cutoff ( Hz )
params.ncomp = 20; % Maximum number of HOSD components
%Params for multivariate HOSD
params.do_multivariate = true; %If the number of channels passed to HOS_spike_detect is greater than 1, use a multivariate version of HOSD, otherwise loop through channels
params.PCA = 0; %Use PCA dimensionality reduction for multivariate data. Values >0 <1 will be treated as proportion of variance to keep, while integers >=1 are taken as the number of PCs to keep.
params.pre_multiply = 1;
params.ncomp_per_chan = 3; % Components per channel for multivariate HOSD. If this value is 0, then keep ncomp fixed at the value in params.ncomp
params.target_srate = 16e3; %( Hz ) If this value is nonempty, greater than zero and less than the original sampling rate, data will be resampled accordingly.
params.get_spike_waveforms = true; %Include all spike waveforms in the output structure. Set to false to avoid prohibitive RAM requirements with large numbers of channels and spikes  
% params.window = 'sasaki'; %Window applied after segmentation of the input
params.window = 'hann'; %Window applied after segmentation of the input
params.tolerance = .001; % Tolerance for the purpose of peak detection (width of smoothing kernel in sec. applied to peak detection output).
params.skewness_threshold = .05; %Discard components whose filter output skewness falls below this threshold. If this is set to 0, then the number
                                 % of retained HOSD components is fixed at ncomp
params.poverlap = .5; %Default window overlap (this may be adjusted down depending on memory limits). Lower values will reduce
                      %memory and cpu requirements, speeding up processing, at the cost of some statistical power.             
params.keepc = 0; %Number of HOSD components to keep. This is set by skewness threshold. 
params.normalize_power = false; %This option to be removed
params.dbt_denoising = true;% Adaptive line noise removal
params.limit_memory = 2e9; %Limit memory usage per component to this value (approximately) in bytes by decreasing the density of estimation windows. Note that for  
                           %multivariate HOSD with many channels memory requirements may be much greater than this.
                           
params.filt_segment = .001; %Duration of filtered data segments to keep. This can be small because the filter concentrates waveform energy in the center of the window
params.hosd_maxiter=25;  %Maximum iterations in HOSD segment alignment
params.online_plot = true; %Show progress in waveform identification during iterative realignment 
params.HOSD_opts = {}; %Put additional arguments for hosobject here as {...,'keyword',[value],...}
% params.HOSD_opts = {'use_adaptive_threshold',false}; %Put additional arguments for hosobject here as {...,'keyword',[value],...}
params.outlier_rejection = @(x)x + 0./((x-nanmean(x(:))).^4./nansum((x(:)-nanmean(x(:))).^4) < .01); %Any sample that accounts for more than 1% of the overall kurtosis should be rejected. 
                                                                                          %This is relatively lenient because we don't want to reject large spikes in recordings with high SNR.
params.check_features = false; %Check whether a component is likely to represent an artifact based on properties of the spectrum.
params.use_components = []; %If non-empty, use only the indicated components for spike detection.

params.output_file.basename = ''; %Filename base (without extension). If nonempty, save outputs to files. A file with the spike detection output is appended with _spike.mat. A file withe the clustering output is appended with _times.mat
params.output_file.basepath = '';
params.output_file.spikes='';
params.output_file.clusters='';

params.subspace_dim = 0;

params.randseed = []; %Random number generator seed 

%%%%
%%%% Below are the default spike sorting parameters contained in the substructure, sorting_params.
%%%%

params.sorting_params.algorithm = @(x,k)mogclust(x,k,'Replicates',25,'Options',struct('MaxIter',500));
params.sorting_params.use_evalclusters = true; 
params.sorting_params.eval_criterion = 'chi2_criterion';  %How to decide on the number of clusters. 'chi2criterion','evalclusters','none'
params.sorting_params.max_cluster_num= 15; %Maximum number of clusters
params.sorting_params.sort_on = 'hosd_filt_segments'; %Sort on concatenated peri-spike-time HOSD filter output. 
% params.sorting_params.sort_on = 'spike_waves'; %Sort on raw spike waveforms 
params.sorting_params.mog_ppweight = .5; %ppweight parameter for chi2_criterion 
params.sorting_params.use_mog_estimates = true; %use_mog_estimates parameter for chi2_criterion
params.sorting_params.make_plots = true; %Plot cluster information
params.sorting_params.isi_qtiles = [.01  .03 .1 .5 .95]; %Get inter-spike intervals at these quantiles
params.sorting_params.cluster_number_penalty = 0.1; %Penalty added to the chi2-KSStat as a fraction of cluster number. 
params.sorting_params.pca_up = 1.5; %Make pca dims used for clustering this multiple of the number of HOSD components

params.plotting.plot_standardized = true; %Plot standardized data for waveforms.
%%%%
%%%% End of default parameter block
%%%%


if nargin >0
    
    params = check_and_assign(paramsin,params);
    
end

function params = check_and_assign(paramsin,params)

%Assign values specified in paramsin to params, otherwise keep defaults. 
%Throw an error if paramsin contains any field not in params.
if isempty(paramsin)
    paramsin = params;
end
fldin = fieldnames(paramsin);    
for k = 1:length(fldin)
    
    switch fldin{k}
        case {'sorting_params','output_file'} %step into structs specified here
            params.(fldin{k}) = check_and_assign( paramsin.(fldin{k}) , params.(fldin{k})); 
        otherwise % lose any default values specified within a substruct. 
            if ~isfield(params,fldin{k})
                error('%s is not a recognized option',fldin{k})      
            else
                params.(fldin{k}) = paramsin.(fldin{k});
            end
    end
end
