% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/misctools/minmax.m $
% $Revision: 193 $
% $Date: 2013-04-28 12:38:09 -0500 (Sun, 28 Apr 2013) $
% $Author: ckovach $
% ------------------------------------------------

function mnmx = minmax(X)

mnmx = cat(1,min(X),max(X));
