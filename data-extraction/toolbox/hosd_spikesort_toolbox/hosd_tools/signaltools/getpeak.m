function ploc = getpeak(f,peaknum);


%finds peak on smooth functions

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

dsdf = diff(sign(diff(f)));

ploc = find(dsdf) + 1;

[sortmagn,sortind] = sort(abs(f(ploc)));

ploc = ploc(flipud(sortind));

if nargin > 1
    
    ploc  = ploc(1:peaknum);
end
    