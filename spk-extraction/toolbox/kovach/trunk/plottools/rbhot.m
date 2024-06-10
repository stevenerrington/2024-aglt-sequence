% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/plottools/rbhot.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

function cmap = rbhot(clip)

if nargin < 1
    clip = 1;
end

ht = hot(64);
ht = ht(1:end-clip,:);

cld = ht(:,[3 2 1]);

cmap = cat(1,flipud(cld),ht);
