function [PEAKXYZs,PKTLMAP] = manualmap(tlXYZs,pkXYZ,label)



% [PEAKXYZs,PKTLMAP] = manualmap(tlXYZs,pkXYZ)
%
% Manually map coordinates in the template tlXYZs to peak coordinates
% pkXYZ. To Link: using shift+mouse select a peak or group of peaks and a group of
% contacts in the template, then right click and select link. To delete
% unwanted or spurious peaks: click delete.

offset = [0 50 -150 0];

global PEAKXYZs

global PKTLMAP

if nargin < 3
    label = '';
end

fig = figure;
% ax = gca;

pkXYZ(:,4) = 1;

mx = max(pkXYZ);
tlXYZs(:,4) = 1;

tlXYZs = tlXYZs+repmat(mx.*(offset~=0)+offset,size(tlXYZs,1),1);


PEAKXYZs = pkXYZ;

PKTLMAP = zeros(size(PEAKXYZs,1),1);

setappdata(fig,'tlXYZ',tlXYZs);
setappdata(fig,'pkXYZ',pkXYZ);

setappdata(fig,'activepk',[]);
setappdata(fig,'activetl',[]);

setappdata(fig,'label',label);
setappdata(fig,'plotstereo',false);

draw(fig);
ax = get(fig,'children');

setappdata(fig,'mainax',ax);
set(ax,'ButtonDownFcn',@axcallback);
setappdata(fig,'plothandles',[]);
 


cxm = uicontextmenu;
setappdata(cxm,'parentax',ax);

m1 = uimenu(cxm,'Label','Delete Selected Peaks','Callback',@deletepeaks);
m2 = uimenu(cxm,'Label','Link','Callback',@linkpeaks);
m3 = uimenu(cxm,'Label','Unlink','Callback',@unlinkpeaks);

set(ax,'UIContextMenu',cxm);

uiwait(fig);




%%%%%%%
function draw(fig)
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/manualmap.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------


global PKTLMAP
global PEAKXYZs

TXY = getappdata(fig,'tlXYZ');
XYZ = PEAKXYZs;

pl = getappdata(fig,'plothandles');
delete(pl(ishandle(pl))), pl = [];

activepk = getappdata(fig,'activepk');
activetl =getappdata(fig,'activetl');
lbl = getappdata(fig,'label');
plst = getappdata(fig,'plotstereo');
if plst
    plotfun = @plotstereo;
else
 plotfun = @plot3;
end
pl = [pl,plotfun(XYZ(:,1),XYZ(:,2),XYZ(:,3),'ko')];

ax = get(pl,'parent');
if iscell(ax),ax=[ax{:}];end

set(ax,'nextplot','add')

% hold off
pl = [pl,plotfun(TXY(:,1),TXY(:,2),TXY(:,3),'+')];

for i = 1:32:size(TXY,1)
    text(TXY(i,1),TXY(i,2),TXY(i,3),sprintf('%i',i))
end


arrayfun(@(ax)grid(ax,'on'),ax)
arrayfun(@(ax)axis(ax,'vis3d'),ax)



% set(pl(1),'zdata',ones(size(TXY,1),1)*(mean(XYZ(:,3))-diff(minmax(XYZ(:,3)))));

% grid on
% axis equal
drawnow

if ~isempty(activepk)
    pl= [pl,plotfun(XYZ(activepk,1),XYZ(activepk,2),XYZ(activepk,3),'r.')];
end
if ~isempty(activetl)
    pl = [pl,plotfun(TXY(activetl,1),TXY(activetl,2),TXY(activetl,3),'r.')];
end

if any(PKTLMAP>0)
    subpk = XYZ(PKTLMAP>0,:);
    subtl = TXY(PKTLMAP(PKTLMAP>0),:);
            
    pl = [pl,plotfun([subpk(:,1),subtl(:,1)]',[subpk(:,2),subtl(:,2)]',[subpk(:,3),subtl(:,3)]','r-')'];

end
setappdata(fig,'plothandles',pl);

set(get(fig,'children'),'ButtonDownFcn',@axcallback);

lbl = getappdata(fig,'label');
title(lbl)

%%%%%%%

function axcallback(ax,e,handles)


fig = get(ax,'parent');

% p1 = get(ax,'currentpoint');
% 
% plhandles = getappdata(fig,'plothandles');
% Xs = get(plhandles,'xdata');
% Ys = get(plhandles,'ydata');
% Zs = get(plhandles,'zdata');
% pkXYZ = cat(2,Xs{1},Ys{1},Zs{1});
% templXYZ = cat(2,Xs{2},Ys{2},Zs{2});
% 
% pkXYZ(:,4) = 1;
% templXYZ(:,4) = 1;
% global PEAKXYZs


switch get(fig,'selectiontype')

%     case 'normal'
      case 'xxx'  %not active, for now
        
       [pkPt, tlPt,mnpk,mntl] = find_near_pt(ax);
       
       if mnpk < mntl
            activept = pkPt;
%             pkactive = true; tlactive = false;
            setappdata(fig,'activepk',activept);
       else
            activetl = tlPt;
%             pkactive = false; tlactive = true;
            setappdata(fig,'activetl',activetl);
       end    
        
        
        
    case 'extend'


%         p2 = get(ax,'currentpoint');
        
       [inpk,intl] = find_within_box(ax);


%         plhandles = getappdata(fig,'plothandles');

        setappdata(fig,'activepk',find(inpk));
        setappdata(fig,'activetl',find(intl));
        draw(fig)

end

    
    
%%%%%% 

function [pkPt, tlPt,mnpk,mntl] = find_near_pt(ax)

global PEAKXYZs

tol = 10;

fig = get(ax,'parent');

p1 = get(ax,'currentpoint');

plhandles = getappdata(fig,'plothandles');
Xs = get(plhandles,'xdata');
Ys = get(plhandles,'ydata');
Zs = get(plhandles,'zdata');


pkXYZ = PEAKXYZs;
templXYZ = getappdata(fig,'tlXYZ');
% templXYZ = cat(1,Xs{2},Ys{2},Zs{2})';
% 
% templXYZ(:,4) = 1;

npks = size(pkXYZ,1);
ncontacts = size(templXYZ,1);

p1(:,4) = 1;


dp1 = diff(p1);
u1 = dp1./sqrt(dp1*dp1');

residmat = eye(4) - u1'*u1;

pkXYZpr = pkXYZ*u1';
pkXYZresid = (pkXYZ - repmat(p1(1,:),npks,1))*residmat;

templXYZpr = templXYZ*u1';
templXYZresid = (templXYZ -repmat(p1(1,:) ,ncontacts,1))*residmat;

dpk = sqrt(sum(pkXYZresid.^2,2));
dtl = sqrt(sum(templXYZresid.^2,2));

[mnpk,pkPt] = min(dpk);
[mntl ,tlPt]= min(dtl);

%%%%%%%%%%

function [inpk,intl] = find_within_box(ax)



global PEAKXYZs

 pkXYZ =PEAKXYZs;
 
fig = get(ax,'parent');
templXYZ = getappdata(fig,'tlXYZ');


p1 = get(ax,'currentpoint');

rect = rbbox;

p2 = get(ax,'currentpoint');

upv = get(ax,'cameraupvector');
upv = upv./sqrt(upv*upv');
upv(4) = 0;


plhandles = getappdata(fig,'plothandles');
% Xs = get(plhandles,'xdata');
% Ys = get(plhandles,'ydata');
% Zs = get(plhandles,'zdata');
% 
% 
% pkXYZ = cat(1,Xs{1},Ys{1},Zs{1});
% templXYZ = cat(1,Xs{2},Ys{2},Zs{2});
% 
% pkXYZ(:,4) = 1;
% templXYZ(:,4) = 1;

npks = size(pkXYZ,1);
ncontacts = size(templXYZ,1);

p1(:,4) = 1;
p2(:,4) = 1;



dp1 = diff(p1);
u1 = dp1./sqrt(dp1*dp1');

p1(1,:) = p1(1,:);
dp12 = (p2(1,:)-p1(1,:))*(eye(4)-u1'*u1);
p2 = p1(1,:)+dp12;

%%% p1 - starting point of box, p2 - ending point, p3 - point at right
%%% angles to the upvector (the lower right corner of the box).
% p3 = p1(1,:) + dp12*(eye(4)-upv'*upv);
p3 = p1(1,:) + dp12*(eye(4)-upv'*upv)*(eye(4)-u1'*u1);
% p4 = p1(1,:) + dp12*upv'*upv;

% u1 = dp1./sqrt(dp1*dp1');


%%% Find transformation from axis space to box-scaled space.
A = [p2(1,:);p3(1,:);p1];
B = [1 1 0 1; 0 1 0 1; 0 0 0 1; 0 0 1 1];

trmat = A\B;

trpk =pkXYZ*trmat;
trtl =templXYZ*trmat;

inpk = all(trpk(:,1:2) > 0 &trpk(:,1:2) < 1,2) ;
intl= all(trtl(:,1:2) > 0 &trtl(:,1:2) < 1,2) ;


%%%%%%%%%
function deletepeaks(h,e,handles)

global PEAKXYZs
global PKTLMAP

ax = getappdata(get(h,'parent'),'parentax');

fig = get(ax,'parent');

% XYZ = getappdata(fig,'pkXYZ');

activepk = getappdata(fig,'activepk');

% XYZ(activepk,:) = [];
PEAKXYZs(activepk,:) = [];
PKTLMAP(activepk) = [];

setappdata(fig,'pkXYZ',PEAKXYZs);
setappdata(fig,'activepk',[]);
% setappdata(ax,'activetl',[]);
set(ax,'nextplot','replace')
draw(fig)

%%%%%%%
function linkpeaks(h,e,handles)

global PKTLMAP
global PEAKXYZs

pkXYZ =PEAKXYZs;

ax = getappdata(get(h,'parent'),'parentax');

fig = get(ax,'parent');

templXYZ = getappdata(fig,'tlXYZ');

activepk = getappdata(fig,'activepk');
% activetl = getappdata(fig,'activetl');

if length(activepk) ~=1 
    setappdata(fig,'activepk',[]);
    return
else

    bdownfun = get(ax,'ButtonDownFcn');
%     bdftemp = @(h,e,hdl) setappdata(h,'dolink',true);
    
    pkxyz = pkXYZ(activepk,:);
    set(gcf,'WindowButtonMotionFcn', @(a,b)plotpos(a,b,pkxyz)); %setting windowButtonMotion Fcn forces update at mouse motion
    
%     setappdata(fig,'dolink',u);
    set(ax,'UserData',[])
    set(gcf,'WindowButtonDownFcn',@(a,b,c) set(ax,'UserData','go'));
    
    waitfor(ax,'UserData','go');
        
    set(ax,'ButtonDownFcn',bdownfun);
    set(gcf,'WindowButtonMotionFcn', []); %setting windowButtonMotion Fcn forces update at mouse motion
    delete(getappdata(fig,'plh'))
    
    [pkPt, tlPt,mnpk,mntl] = find_near_pt(ax);
    tol = 10;
    
    if mntl < tol
        PKTLMAP(activepk) = tlPt;
    end
    
    draw(fig)
    
    
end

%%%%%%%%%%%%%%%%%

function plotpos(fig,e,pkxyz)

ax = getappdata(fig,'mainax');

fig = get(ax,'parent');

cp = get(ax,'currentpoint');

plh = getappdata(fig,'plh');

if isempty(plh) || ~ishandle(plh)
    plh = plot3([pkxyz(1),cp(1,1)],[pkxyz(2),cp(1,2)],[pkxyz(3),cp(1,3)],'r--');
    setappdata(fig,'plh',plh)
else
    set(plh,'xdata',[pkxyz(1),cp(1,1)],'ydata',[pkxyz(2),cp(1,2)],'zdata',[pkxyz(3),cp(1,3)]);
end


%%%%%%%
function unlinkpeaks(h,e,handles)

global PKTLMAP

ax = getappdata(get(h,'parent'),'parentax');

fig = get(ax,'parent');

activepk = getappdata(fig,'activepk');
activetl = getappdata(fig,'activetl');

if length(activepk) ~=1 
    setappdata(fig,'activepk',[]);
    setappdata(fig,'activetl',[]);
    return
else

    PKTLMAP(activepk) = 0;
end

