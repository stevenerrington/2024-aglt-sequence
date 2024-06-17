
function [cl,gm,gms]=mogclust(x,n,varargin)


gms = {};
for k = 1:length(n)
    gms{k} = fitgmdist(x,n(k),varargin{:});
end

BIC = cellfun(@(x)x.BIC,gms);
gm=gms{BIC==min(BIC)};
cl = gm.cluster(x);
