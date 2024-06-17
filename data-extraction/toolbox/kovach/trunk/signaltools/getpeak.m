function ploc = getpeak(f,peaknum);


%finds peak on smooth functions

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/signaltools/getpeak.m $
% $Revision: 676 $
% $Date: 2016-01-15 14:07:43 -0600 (Fri, 15 Jan 2016) $
% $Author: ckovach $
% ------------------------------------------------

dsdf = diff(sign(diff(f)));

ploc = find(dsdf) + 1;

[sortmagn,sortind] = sort(abs(f(ploc)));

ploc = ploc(flipud(sortind));

if nargin > 1
    
    ploc  = ploc(1:peaknum);
end
    