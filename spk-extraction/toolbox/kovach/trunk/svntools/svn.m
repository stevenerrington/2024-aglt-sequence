function varargout = svn(varargin)

%
% Pass svn commands to the system. This is equivalent to !svn ...
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/svntools/svn.m $
% $Revision: 285 $
% $Date: 2013-07-29 08:33:45 -0500 (Mon, 29 Jul 2013) $
% $Author: ckovach $
% ------------------------------------------------

com = ['svn ',sprintf('%s ',varargin{:})];

[stat,res] = system(com);
fprintf(res)
if nargout > 0
    varargout{1} = stat;
end
if nargout > 1 
    varargout{2} = res;
end

