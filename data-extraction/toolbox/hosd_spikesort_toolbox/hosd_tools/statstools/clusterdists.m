
function [c2c_dists, dists2center,clust_id,project,discarded,SS,M,Chi2_KSstat] = clusterdists(x,c,threshold,M,SS)

%%% Returns squared Mahalanobis distances between cluster centers and
%%% within cluster distance to center for each sample.


if nargin < 3
    threshold = 5;
end
[clust_id,~,unqi] = unique(c);

project = [];

dists2center = nan(size(x,1),1);
c2c_dists=[];

nclust = length(clust_id);
sqf = squareform(1:nclust*(nclust-1)/2);
for k = find(~isnan(clust_id(:)))'
   
    run = true;
    while run
        geti = find(unqi==k);
        m = mean(x(geti,:));
        xc = x(geti,:)-m;
        N = size(xc,1);
        
        if nargin < 5
            SS(:,:,k) = xc'*xc /(N-1);  
        end
        if nargin < 4
            M(k,:) = m;
        end
        
        d2c = sum((xc*SS(:,:,k)^-1).*xc,2);%*sum(unqi==k);
        discard = d2c>threshold.^2;
        run = any(discard);
        unqi(geti(discard))=0;
    end
    
    dists2center(unqi==k,:) = d2c;
    
    for kk = 1:k-1
       
        mS = (SS(:,:,k) + SS(:,:,kk));
        dm = M(k,:)-M(kk,:);
        
%         c2c_dists(end+1,:) = dm*mS^-1*dm';
        c2c_dists(sqf(k,kk),:) = dm*mS^-1*dm'; %#ok<*AGROW>
        
        project(:,sqf(k,kk)) = mS^-.5*dm'./sqrt(sum(dm.^2));
    end
end



if nargout> 7
    %For Gaussian clusters, this should be uniformly distributed
    u = chi2cdf(dists2center,size(x,2));
    Chi2_KSstat = max(abs(sort(u(:))-linspace(0,1,length(u))'));
end

discarded = unqi==0;
