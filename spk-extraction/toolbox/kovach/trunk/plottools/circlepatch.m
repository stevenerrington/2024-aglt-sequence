
function pt = circlepatch(r,cx,cy,c,alpha,ax,scale)

%circlepatch(r,cx,cy,c)
% Draws circular patches on the axis, where r is radius, cx is center xs,
% cy is center ys, and c is index into colormap for eatch patch

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/plottools/circlepatch.m $
% $Revision: 724 $
% $Date: 2016-04-07 20:50:13 -0500 (Thu, 07 Apr 2016) $
% $Author: ckovach $
% ------------------------------------------------

nsides =20;
if nargin < 5
    alpha = 1;
end
if nargin< 6
    ax = gca;
end
if nargin< 7
    scale = [1 1];
end

if nargin < 3 || isempty(cx) || isempty(cy)
    axes(ax)
    [cx,cy] = ginput(1);
end
    
if length(r) == 1
    r = r*ones(size(cx));
end

axes(ax)

dth = 2*pi./nsides;
th = (dth:dth:2*pi)';

x = scale(1)*cos(th)*r(:)' + ones(size(th))*cx(:)';
y = scale(2)*sin(th)*r(:)' + ones(size(th))*cy(:)';


pt = patch(x,y,repmat(c,length(th),1));

set(pt,'facealpha',alpha)
