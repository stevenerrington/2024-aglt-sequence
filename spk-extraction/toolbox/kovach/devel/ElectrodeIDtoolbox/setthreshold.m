
function sdthresh = setthreshold(pks,pkval)

global thrval

slidemax = 20;
slidemin = 0;

fig = figure;


ui = uicontrol('style','slider',...
                'Units','normalized',...
                'position',[0 0 1 .1],...
                'min',slidemin,'max',slidemax,'value',2);
            

thr = get(ui,'value');
plpk = pks; plpk(pkval< thr,:) = nan;            
pl = plot3(plpk(:,1),plpk(:,2),plpk(:,3),'.');
axis equal
grid on
title(sprintf('%f',thr));

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/setthreshold.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------


setappdata(ui,'pks',pks)
setappdata(ui,'pkval',pkval)
setappdata(ui,'pl',pl)
setappdata(ui,'parent',fig)
set(ui,'Callback',@ uicallback)
uiwait(fig)

sdthresh =thrval;

%%%%


function uicallback(ui,evt,h)

global thrval

pks = getappdata(ui,'pks');
pl = getappdata(ui,'pl');
pkval = getappdata(ui,'pkval');
fig = getappdata(ui,'parent');

thrval = get(ui,'value');
title(sprintf('%f',thrval));

plpk = pks; plpk(pkval< thrval,:) = nan;         

set(pl,'Xdata',plpk(:,1)','Ydata',plpk(:,2)','Zdata',plpk(:,3)')
drawnow


