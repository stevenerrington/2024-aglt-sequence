
function [rgb,ind] = cm2rgb(val,fig,ax)

%returns the rgb value and index for a value using the current colormap.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 3
    fig = gcf;
    ax = gca;
end

if ~isscalar(fig) || ~ishandle(fig) && isnumeric(fig)
    cmap = fig;
    cax = ax;
else    
    cmap = get(fig,'colormap');
    cax = caxis(ax);
end

ind = ceil( (val-cax(1))./range(cax)*size(cmap,1));

ind(ind<=0|isnan(ind)) = 1;
ind(ind>size(cmap,1)) = size(cmap,1);

rgb = cmap(ind(:),:);


