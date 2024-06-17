function clout = clusterfork(x,threshold,clusterfun,clin)

% clout = clusterfork(x,threshold,clusterfun,clin)
%
% Each cluster forks itself, yielding two offspring, until the Mahalanobis distance bewtween offspring
% falls below a threshold.
%
% Input: 
%   x - input data with rows as observations
%   threshold - Mahalanobis threshold (default = 4);
%   clusterfun - cluster function, which must take 2 arguments: data and
%                number of clusters. Defaut:  @(x,k)kmeans(x,k,'maxiter',500,'replicates',20);
%
% C. Kovach 2020

nchild = 2; % Number of offspring produced by each fork.

if nargin < 2 || isempty(threshold)
    threshold = 4;
end

if nargin < 3 || isempty(clusterfun)
    clusterfun = @(x,k)kmeans(x,k,'maxiter',500,'replicates',20);
end

if nargin < 4 || isempty(clin)
    clin = zeros(size(x,1),1);
end

minsize = size(x,2);

clout = clin;



[unq,~,unqi] = unique(clin);

for k = 1:length(unq)
    
    xx = x(unqi==k,:);
    

   cl = clusterfun(xx,nchild);
   [~,~,cl] = unique(cl);
   cl = cl-1;
   
    if min(crosstab(cl))<minsize
     continue
    end   
    
   cldist = clusterdists(xx,cl);
   if ~isnan(cldist) && sqrt(cldist)>=threshold        
       cl = cl + nchild*(clusterfork(xx,threshold,clusterfun,cl)-1);
       clout(unqi==k) = cl;
   end
   
end

[~,~,clout] = unique(clout);

