function varargout=reseed(s)

%
% RESEED seeds the uniform and normal random generators, RAND and RANDN
% with the current time obtained with NOW.
%
% See also: RAND, RANDN, NOW.

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/statstools/reseed.m $
% $Revision: 1140 $
% $Date: 2019-01-08 11:35:40 -0600 (Tue, 08 Jan 2019) $
% $Author: ckovach $
% ------------------------------------------------

umethod = 'twister';
nmethod = 'state';

if nargin < 1
    s = typecast(now,'uint32');
    s=s(1); %Using the least significant 4 bytes of the value returned by now
end

rand(umethod, double( s )); 

srandn=double( swapbytes( s ) );
randn( nmethod, srandn ); %using swapbytes so that the randn and rand seeds are
                                             %at least somewhat arbitrary
                                             %with respect to each other.
                                             
if nargout > 0
    varargout{1}=s;
end

