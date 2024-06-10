function rotor(indx)

%Displays a rotating icon inside a loop
r = '|/-\';
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/misctools/rotor.m $
% $Revision: 36 $
% $Date: 2011-06-04 23:10:51 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------


if indx > 1
    fprintf('\b')
end

fprintf('%s',r(mod(indx,4)+1));

