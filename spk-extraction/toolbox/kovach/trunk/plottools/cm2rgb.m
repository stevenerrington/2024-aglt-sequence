
function [rgb,ind] = cm2rgb(val,fig,ax)

%returns the rgb value and index for a value using the current colormap.

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/plottools/cm2rgb.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 3
    fig = gcf;
    ax = gca;
end

cmap = get(fig,'colormap');
cax = caxis(ax);

ind = ceil( (val-cax(1))./range(cax)*size(cmap,1));

ind(ind<=0) = 1;
ind(ind>size(cmap,1)) = size(cmap,1);

rgb = cmap(ind(:),:);


