function setstereo(ax1,ax2,sep)

% setstero(ax1,ax2,[sep]). Locks the positions of the cameras in two axes to
% provide a stereoscopic view. Interocular distance is given in 'sep' as a
% proportion of the distance to the target (default, sep = .05).
%
% See also UPDATESTEREO and PLOTSTEREO

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/plottools/setstereo.m $
% $Revision: 935 $
% $Date: 2017-10-26 15:06:53 -0500 (Thu, 26 Oct 2017) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 3
    sep = .05;
end

fg1 = get(ax1,'parent');
fg2 = get(ax2,'parent');

set(ax1,'PlotBoxAspectRatio',[1 1 1])
set(ax2,'PlotBoxAspectRatio',[1 1 1])
set(ax2,'cameraviewangle',get(ax1,'cameraviewangle'))
set(ax2,'dataaspectratio',get(ax1,'dataaspectratio'))
set([ax1,ax2],'dataaspectratiomode','manual','plotboxaspectratiomode','manual','cameraviewanglemode','manual')
updatestereo(ax1,ax2,sep);
set(fg1,'WindowButtonMotionFcn',@(a,b)updatestereo(ax1,ax2,sep) )
set(fg2,'WindowButtonMotionFcn',@(a,b)updatestereo(ax1,ax2,sep) )
