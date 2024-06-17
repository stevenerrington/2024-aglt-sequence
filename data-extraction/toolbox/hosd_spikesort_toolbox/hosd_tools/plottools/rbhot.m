% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

function cmap = rbhot(clip,n)

if nargin < 1 || isempty(clip)
    clip = 1;
end
if nargin < 2 || isempty(n)
    n=64;
end
ht = hot(floor((n+clip)/2));
ht = ht(1:end-floor(clip/2),:);

cld = ht(:,[3 2 1]);

cmap = cat(1,flipud(cld),ht);
