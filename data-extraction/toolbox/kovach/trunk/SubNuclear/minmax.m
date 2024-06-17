% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/SubNuclear/minmax.m $
% $Revision: 1152 $
% $Date: 2019-02-03 21:24:39 -0600 (Sun, 03 Feb 2019) $
% $Author: ckovach $
% ------------------------------------------------

function [mnmx,mnmxi] = minmax(X)

if nargout <2
    mnmx = cat(1,min(X),max(X));
else
   [mn,mni] = min(X);
   [mx,mxi] = max(X);
   mnmx = [mn;mx];
   mnmxi = [mni;mxi];
end