
function [jkm,jkse,cohjk, csjk,sp1jk,sp2jk] = cohjacknife(f1,f2,ewgts)

%Computes jacknife estimate of variance for the coherence and crosspectra
%of two complex time-frequency decomposed signals f1 and f2

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/cohjacknife.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

if size(f1,3)>1 || size(f2,3)>1
    dim = 3;
else
    dim = 2;
end

if nargin < 3
    ewgts = 1;
end

n = size(f1,dim);

rep = ones(1,3); rep(dim) = n;

repfun = @(x,y) repmat(sum((x.*conj(y))*ewgts,dim),rep) - x.*conj(y);

csjk = repfun(f1,f2);
sp1jk = repfun(f1,f1);
sp2jk = repfun(f2,f2);

cohjk = csjk./sqrt(sp1jk.*sp2jk);

jkm = mean(abs(cohjk),dim);

jkse = sqrt((mean(abs(cohjk).^2,dim) - jkm.^2)*(n-1));


