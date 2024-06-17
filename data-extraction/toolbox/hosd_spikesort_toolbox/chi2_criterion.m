
function [out,cleval,CL] = chi2_criterion(x,nclusts,clusterfun,use_mog_estimates ,ppweight,min_cluster_size,pool_dists,cluster_number_penalty)

%  [out,cleval] = chi2_criterion(X,nclusts,[clusterfun],[use_mog_estimates] ,[ppweight])
%
% Choose the optimal cluster number according to the deviation of cluster distances
% from a chi-square distribution. The rationale is to find clustering 
% that creates the most Gaussian-like within-cluster distribution.
% Kolmogorov-Smirnov statistic is used as the criterion value.
%
% Inputs:
%   X - p-dimensional data as an N x P matrix
%   nclusts - Cluster numbers to evaluate. If nclusts is scalar, then
%             1:nclusts will be evaluated.
%   clusterfun - clustering function. Default is a mixture of Gaussians
%                 with 10 replicates.
%   use_mog_estimates - If the second output of clusterfun is a gmdistribution object,
%                 then the mixture of Gaussians estimates for cluster means
%                 and variances will be used rather than sample
%                 means if use_mog_estimates is true. Otherwise
%                 sample values are used. Default is true.
%   ppweight -  When clustering is based on a mixture of Gaussians, this
%               parameter selects between three possibilities for how 
%               distances can be assigned to samples:
%                  1) According to cluster assignment (ppweight = 0)
%                  2) Randomly according to the posterior probability of 
%                     source distribution (ppweight = 1).
%                  3) Compromise between the two: instead of posterior probability,
%                     assign randomly according to
%                       P = ppweight*posterior_probability + (1-ppweight)*cluster_assignment
%               In case 1, the criterion tends to over-penalize cluster
%                   number causing overlapping clusters to be lumped together.
%               In case 2, the criterion tends to under-penalize cluster
%                   number leading to inappropriate splitting of clusters.
%               Case 3, provides a way to find a suitable middle ground.
%               Default value is ppweight = .75.
%

% C. Kovach 2021

look_ahead = 3; %Stop if minimum value hasn't been attained within the most recent look_ahead steps. 

if nargin < 3 || isempty(clusterfun)
    clusterfun = @(x,k)mogclust(x,k,'Replicates',10,'Options',struct('MaxIter',500));
end

if nargin < 4 || isempty(use_mog_estimates )
    use_mog_estimates = true; % Using the mean and variance estimates from the gm fitting tends to
end                        % generate fewer clusters.

if nargin < 5 || isempty(ppweight)
    ppweight = .75; % Compromise between posterior probability and cluster assignment in randomly assigning distance to points 
                    % 1 tends to under penalize cluster number while 0 tends to over penalize.
end

if ~use_mog_estimates
    ppweight = nan;
end

if isscalar(nclusts)
    nclusts = 1:nclusts;
end
if nargin < 6 || isempty(min_cluster_size)
    min_cluster_size = 5; %Reject if clustering produces one or more clusters with fewer members than this value.
end

if nargin < 7 || isempty(pool_dists)
    pool_dists = false; % If true compute KS stats after pooling distances from all clusters. Otherwise compute per cluster and average the result.
end
if nargin < 8 || isempty(cluster_number_penalty)
    cluster_number_penalty =  0; %An additional penalty for cluster number
end

fprintf('\nEvaluating clusters using the chi-square criterion.')
for k = 1:length(nclusts)
    
    fprintf('\n\nTrying %i clusters with %s...',nclusts(k),char(clusterfun))
           
    
    tic
    try
        if contains(char(clusterfun),'mogclust')
            [cl,gm] = clusterfun(x,nclusts(k));
            fprintf('\nCompleted in %f sec',toc)
            cleval(k).gmdist = gm;
            cleval(k).gmAIC = gm.AIC;
            cleval(k).gmBIC = gm.BIC;
        else
            cl = clusterfun(x,nclusts(k));
            cleval(k).gmdist.posterior = @(x)nan;
%             cleval(k).gmdist.AIC = nan;
%             cleval(k).gmdist.BIC = nan;
%             cleval(k).gmAIC = nan;
%             cleval(k).gmBIC = nan;
%             cleval(k).gmdist.cluster = @(x)nan(size(x,1),1);
%             cleval(k).gmdist.mahal = @(x)nan(size(x,1),1);
        end
         CL(:,k) = cl;
       
    catch err
        fprintf('\nClustering FAILED with error: %s',err.message)
%         KSstat(k) = nan;
         CL(:,k) = nan;
        cleval(k).gmdist.posterior = @(x)nan;
        cleval(k).gmdist.AIC = nan;
        cleval(k).gmdist.BIC = nan;
        cleval(k).gmAIC = nan;
        cleval(k).gmBIC = nan;
        cleval(k).gmdist.cluster = @(x)nan(size(x,1),1);
        cleval(k).gmdist.mahal = @(x)nan(size(x,1),1);
        cl = [];
%         U(:,k,:,:) = nan;
%         continue
    end      
    for kk = 1:length(ppweight)
        if ~isnan(ppweight(kk)) && any(contains(fieldnames(cleval(k).gmdist),'BIC') )
            mhd = cleval(k).gmdist.mahal(x);
            clustid = cleval(k).gmdist.cluster(x) == (1:nclusts(k));
                post = [zeros(size(x,1),1), cleval(k).gmdist.posterior(x)*ppweight(kk) + (1-ppweight(kk))*clustid];
           
    %          post = chi2pdf(mhd,size(x,2));
    %          post = [zeros(size(x,1),1),post./sum(post,2)];
              r = rand(size(post,1),1);
              rpick = diff(r<cumsum(post,2),[],2);
            cleval(k).d2center(:,kk) = sum(rpick.*mhd,2); %#ok<*SAGROW>
        else
            args = {};
            [cleval(k).c2c,cleval(k).d2center(:,kk)] = clusterdists(x,cl,Inf,args{:}); %#ok<*AGROW>           
            
        end
    end
    [unq,~,unqi] = unique(cl);

    %For Gaussian clusters, this should be uniformly distributed
    if pool_dists
        u = chi2cdf(cleval(k).d2center,size(x,2));
        KSstat(k,:) = max(abs(sort(u)-linspace(0,1,size(u,1))'))*sqrt(size(u,1)) + cluster_number_penalty*nclusts(k);
%         U(:,k,:,:) = u;
    else
        ksst = nan(size(unq));
        for kk = 1:length(unq)
           u =  chi2cdf(cleval(k).d2center(unqi==kk,:),size(x,2));
           ksst(kk,:) = max(abs(sort(u)-linspace(0,1,size(u,1))'))*sqrt(sum(unqi==kk));
        end
        KSstat(k,:) = nanmean(ksst,1) + cluster_number_penalty*nclusts(k); 
    end
    [~,mni] = min(KSstat);
    fprintf('\nChi-square KS stat (ppweight: %0.2f): %0.2f',ppweight(1),KSstat(k,1))
    if k-max(mni) > look_ahead
        fprintf('\nKSstat has not attained a minimum value in %i steps. Stopping.',look_ahead);
        break
    end
    if length(unq)< nclusts(k) || any(hist(unqi,1:length(unq))<min_cluster_size)
        fprintf('\nExcluding because one or more clusters fell below min. size of %i',min_cluster_size) 
        KSstat(k,:) = Inf;
    end
end


[~,out.OptimalK] = min(KSstat); 
out.OptimalY = CL((1:size(CL,1))' + (out.OptimalK-1)*size(CL,1));
out.KSstat = KSstat;
if isfield(cleval,'gmBIC')
    out.BIC = [cleval.gmBIC];
end
out.nclusts= nclusts;
out.eval = cleval(out.OptimalK);
out.ppweight = ppweight;
