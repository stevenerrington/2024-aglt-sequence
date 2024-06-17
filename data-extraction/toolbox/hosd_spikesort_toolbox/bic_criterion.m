
function [out,cleval,CL] = bic_criterion(x,nclusts,clusterfun,min_cluster_size,pool_dists,cluster_number_penalty)

%  [out,cleval] = bic_criterion(X,nclusts,[clusterfun],[use_mog_estimates] ,[ppweight])
%
% Choose the optimal cluster number for a mixture of Gaussians according to BIC. 
%
% Inputs:
%   X - p-dimensional data as an N x P matrix
%   nclusts - Cluster numbers to evaluate. If nclusts is scalar, then
%             1:nclusts will be evaluated.
%   clusterfun - clustering function. Default is a mixture of Gaussians
%                 with 10 replicates.


% C. Kovach 2021

look_ahead = 5; %Stop if minimum value hasn't been attained within the most recent look_ahead steps. 

if nargin < 3 || isempty(clusterfun)
    clusterfun = @(x,k)mogclust(x,k,'Replicates',10,'Options',struct('MaxIter',500));
end


if isscalar(nclusts)
    nclusts = 1:nclusts;
end

if nargin < 4 || isempty(min_cluster_size)
    min_cluster_size = 5; %Reject if clustering produces one or more clusters with fewer members than this value.
end

if nargin < 5 || isempty(pool_dists)
    pool_dists = false; % If true compute KS stats after pooling distances from all clusters. Otherwise compute per cluster and average the result.
end
% if nargin < 6 || isempty(cluster_number_penalty)
%     cluster_number_penalty =  0; %An additional penalty for cluster number
% end

fprintf('\nEvaluating clusters using BIC.')
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
            error('Must use mogclust')
        end
         CL(:,k) = cl;
       
    catch err
        fprintf('\nClustering FAILED with error: %s',err.message)
%         BIC(k) = nan;
         CL(:,k) = nan;
        cleval(k).gmdist.posterior = @(x)nan;
        cleval(k).gmdist.AIC = nan;
        cleval(k).gmdist.BIC = nan;
        cleval(k).gmAIC = nan;
        cleval(k).gmBIC = nan;
        cleval(k).gmdist.cluster = @(x)nan(size(x,1),1);
        cleval(k).gmdist.mahal = @(x)nan(size(x,1),1);
%         U(:,k,:,:) = nan;
%         continue
    end      
%     for kk = 1:length(ppweight)
%         if ~isnan(ppweight(kk)) && any(contains(fieldnames(cleval(k).gmdist),'BIC') )
%             mhd = cleval(k).gmdist.mahal(x);
%             clustid = cleval(k).gmdist.cluster(x) == (1:nclusts(k));
%                 post = [zeros(size(x,1),1), cleval(k).gmdist.posterior(x)*ppweight(kk) + (1-ppweight(kk))*clustid];
%            
%     %          post = chi2pdf(mhd,size(x,2));
%     %          post = [zeros(size(x,1),1),post./sum(post,2)];
%               r = rand(size(post,1),1);
%               rpick = diff(r<cumsum(post,2),[],2);
%             cleval(k).d2center(:,kk) = sum(rpick.*mhd,2); %#ok<*SAGROW>
%         else
%             args = {};
%             [cleval(k).c2c,cleval(k).d2center(:,kk)] = clusterdists(x,cl,Inf,args{:}); %#ok<*AGROW>           
%             
%         end
%     end
    [unq,~,unqi] = unique(cl);

    %For Gaussian clusters, this should be uniformly distributed
    BIC = [cleval.gmBIC];
    [~,mni] = min(BIC);
    fprintf('\nChi-square Delta-BIC: %0.2f',BIC(k)-min(BIC(1:k-1)))
    if k-mni > look_ahead
        fprintf('\nBIC has not attained a minimum value in %i steps. Stopping.',look_ahead);
        break
    end
    if length(unq)< nclusts(k) || any(hist(unqi,1:length(unq))<min_cluster_size)
        fprintf('\nExcluding because one or more clusters fell below min. size of %i',min_cluster_size) 
        BIC(k) = Inf;
    end
end


[~,out.OptimalK] = min(BIC); 
out.OptimalY = CL((1:size(CL,1))' + (out.OptimalK-1)*size(CL,1));
out.BIC = BIC;
if isfield(cleval,'gmBIC')
    out.BIC = [cleval.gmBIC];
end
out.nclusts= nclusts;
out.eval = cleval(out.OptimalK);
