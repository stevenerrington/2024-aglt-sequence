function [P,T,df,Q] = polyt(X,Y)

% [P,Tsq] = polyt(X,Y)
% Perfoms multiple uncorrected t-tests on the rows of X an Y, which must be the same
% length.

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/statstools/polyt.m $
% $Revision: 56 $
% $Date: 2011-09-28 12:41:31 -0500 (Wed, 28 Sep 2011) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 2
    Y = 0;
end
dim = 2;
if size(X,3)>1 || size(Y,3)>1
    dim = 3;
end

if prod(size(Y)) == 1;
    onegr = 1;
else
   onegr = 0; 
end


nx = size(X,dim);
ny = size(Y,dim);

mx = mean(X,dim);
my = mean(Y,dim);

if dim < 3
    SSx = sum((X- mx*ones(1,size(X,2))).^2,dim);
    SSy = sum((Y- my*ones(1,size(Y,2))).^2,dim);
else
    SSx = sum((X- repmat(mx,[1 1 nx])).^2,dim);
    SSy = sum((Y- repmat(my,[1 1 ny])).^2,dim);
end    

if onegr 
    serm = sqrt(SSx./(nx*(nx-1)));
    T = (mx - Y)./serm;
    df = nx - 1;
else
    Spl = (SSx+SSy)/(nx+ny-2);

    Tsq = nx*ny/(nx+ny)*((mx-my).^2)./Spl;

    T = sqrt(Tsq).*sign(mx-my);
    df = (nx+ny-2)*ones(size(mx));
end
p = tcdf(T,df);


P = 2*min(cat(dim,p,1-p),[],dim);  %2-tialed p;

if nargout > 3
    Q = nan(size(P));
   Q(:) = fdr(P(:)) ;
end