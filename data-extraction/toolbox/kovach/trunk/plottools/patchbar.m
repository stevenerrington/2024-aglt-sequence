function h = patchbar(x,y,cdata,width,baseline)


% patchbar(x,y,cdata,[width])
%
% Generates a bar plot with patch so that the face color of each bar can be
% be set separately.
%
% Input:
%
% x - bar x position
% y - bar height
% cdata - color data for each bar
% width - bar width. Defaults to 2/3 avg spacing of x.

%
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/plottools/patchbar.m $
% $Revision: 934 $
% $Date: 2017-10-26 15:06:28 -0500 (Thu, 26 Oct 2017) $
% $Author: ckovach $
% ------------------------------------------------
%
% C. Kovach 2015


if nargin < 3 || isempty(cdata)
    cdata = 'b';
end

if isempty(x)   
    x = 1:length(y);
end
if nargin < 4 || isempty(width)
    width = 2/3*abs(mean(diff(x)));
end
if nargin < 5 || isempty(baseline)
    baseline=0;
end

X = [x-width/2;x+width/2;x+width/2;x-width/2];
Y = [0 0 1 1]'*y(:)';
if baseline~=0
    Y = Y+baseline;
end

h = patch(X,Y,[1 1 1 1]'*cdata(:)');
