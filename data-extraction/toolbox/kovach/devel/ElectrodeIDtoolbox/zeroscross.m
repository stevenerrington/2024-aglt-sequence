function PK = zeroscross(X,sign)

%find the intersection of zero crossings of diff(X) along all dimensions
%(i.e. find the peaks).

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/zeroscross.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------


if nargin < 2
    sign = 1; %Find peaks if sign  ==1 , valleys if sign == -1;
end
sz = size(X);
ndim = length(sz);

if ndim == 2 && sz(2) == 1, ndim = 1; end

Dfun = @(Z,i) diff(Z,[],i);
% D2 = @(Z,i,j) D(D(Z,[],i),j);

indxfun = @(sti,endi) cellfun(@(a,b) a:b,num2cell(sti),num2cell(endi),'uniformoutput',false);


ZC = zeros(sz);

for i = 1:ndim
    
    
    z = zeros(1,ndim);
    z(i) = 1;

    d = Dfun(X,i);
    
    a = indxfun(1+z,sz-z);
    b = indxfun(1+z*0,sz-z*2);
    c = indxfun(1+z,sz-z);
    ZC(a{:}) = ZC(a{:}) + ((sign*d(b{:}) > 0) & (sign*d(c{:}) < 0));  
end
PK = ZC==ndim; %Points which are zero crossings along each of the dimensions are identified as peaks


