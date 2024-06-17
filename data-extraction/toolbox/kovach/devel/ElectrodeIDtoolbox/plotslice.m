function plotslice(IM,pts,hdr,decim)


% function plotslice(IM,pts,hdr)
% 
%  A function to plot x-y-z- slices of imaging data overlayed with point data. 
%
%   IM - Image
%   pts - points to plot in image
%   hdr - header data as returned by readnifit. 
%        OR transformation matrix to be applied to voxel units (default is eye(4)).
%
% C. Kovach 2010
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/plotslice.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

global IMNRM ZDEC XDEC YDEC



if nargin< 3
    tmat =eye(4);
end

if nargin < 4 || isempty(decim)
    decim = 1;
end

imdim = size(IM);

[I,J,K] = meshgrid(1:imdim(2),1:imdim(1),1:imdim(3));
IJK = [I(:),J(:),K(:),ones(numel(I),1)];
IJK(:,4) = 1;

if nargin < 3
    trm = eye(4);
elseif isstruct(hdr)
    trm = double(hdr.vox2unit');
    trm(4,4) = 1;
else
    trm = hdr;
end



XYZ = IJK*trm;
mn = min(XYZ);
mx = max(XYZ);



X = reshape(XYZ(:,1),size(I));
% x = squeeze(X(1,1:decim:end,1));
Y = reshape(XYZ(:,2),size(I));
% y = squeeze(Y(1:decim:end,1,1));
Z = reshape(XYZ(:,3),size(I));
% z = squeeze(Z(1,1,1:decim:end));
ZDEC = Z(1:decim:end,1:decim:end,1:decim:end);
XDEC = X(1:decim:end,1:decim:end,1:decim:end);
YDEC = Y(1:decim:end,1:decim:end,1:decim:end);
IMdecim = IM(1:decim:end,1:decim:end,1:decim:end);
IMNRM = (IMdecim-min(IMdecim(:)))./range(IMdecim(:));


fig = figure;
ax = axes;

set(fig,'units','normalized','DeleteFcn',@cleanup)
set(ax,'units','normalized','position',[.02 .02 .75 .96]);
setappdata(fig,'imhandle',nan(1,3));


setappdata(fig,'IMNRM',IMNRM);
setappdata(fig,'XDEC',XDEC);
setappdata(fig,'YDEC',YDEC);
setappdata(fig,'ZDEC',ZDEC);

axlim(1:2:5) = mn(1:3);
axlim(2:2:6) = mx(1:3);

setappdata(fig,'axlim',axlim);


% setappdata(fig,'sl',);
% setappdata(fig,'sl',ceil(imdim(2)/2));
% setappdata(fig,'sl',ceil(imdim(3)/2));

% Create user interface (check box and sliders)

hd.xcheck = uicontrol(fig,'Style','checkbox','units','normalized','position',[.6 .95 .1 .05],'String','X: .5','Callback',@redrawx,'Value',1);
hd.ycheck = uicontrol(fig,'Style','checkbox','units','normalized','position',[.6 .9 .1 .05],'String','Y: .5','Callback',@redrawy,'Value',1);
hd.zcheck = uicontrol(fig,'Style','checkbox','units','normalized','position',[.6 .85 .1 .05],'String','Z: .5','Callback',@redrawz,'Value',1);
hd.xslide = uicontrol(fig,'Style','slider','units','normalized','position',[.7 .95 .3 .05],'Callback',@redrawx,'Value',.5);
hd.yslide = uicontrol(fig,'Style','slider','units','normalized','position',[.7 .9 .3 .05],'Callback',@redrawy,'Value',.5);
hd.zslide = uicontrol(fig,'Style','slider','units','normalized','position',[.7 .85 .3 .05],'Callback',@redrawz,'Value',.5);
hd.aslide = uicontrol(fig,'Style','slider','units','normalized','position',[.7 .8 .3 .05],'Callback',@redraw,'Value',.5);
hd.alphalbl = uicontrol(fig,'Style','text','units','normalized','position',[.58 .8 .1 .05],'String','alpha: .5','Value',.5);

if ~isempty(pts)
    pkpl = plot3(pts(:,1),pts(:,2),pts(:,3),'r.');
else
    pkpl = [];
end
setappdata(fig,'pkh',pkpl);
hold on,
axis(axlim);
axis equal
grid on



setappdata(fig,'handles',hd);

redraw(ax);

view(3);

%%%%%%%%%%%%%

function redrawx(h,ev,hd)

%Redraw X slice

global IMNRM ZDEC XDEC YDEC

imdim = size(IMNRM);
fig = get(h,'parent');

hd = getappdata(fig,'handles');

imh = getappdata(fig,'imhandle');


if get(hd.xcheck,'Value')
    sl = ceil(get(hd.xslide,'Value')*imdim(1));   
      if ~ishandle(imh(1))
           imh(1) =surf(squeeze(XDEC(sl,:,:)),squeeze(YDEC(sl,:,:)),squeeze(ZDEC(sl,:,:)),squeeze(IMNRM(sl,:,:)));
    else
        set(imh(1),'Xdata',squeeze(XDEC(sl,:,:)),'Ydata',squeeze(YDEC(sl,:,:)),'Zdata',squeeze(ZDEC(sl,:,:)),'Cdata',squeeze(IMNRM(sl,:,:)));
    end
      set(imh(1),'facealpha',get(hd.aslide,'value'));    
      set(hd.xcheck,'String',sprintf('X: %3i',sl))
elseif ishandle(imh(1)) 
    delete(imh(1));
          set(hd.xcheck,'String',sprintf('X: %3i',[]))

end
setappdata(fig,'imhandle',imh);
shading flat


%%%%%%%%%%%%%

function redrawy(h,ev,hd)
%Redraw Y slice


global IMNRM ZDEC XDEC YDEC

imdim = size(IMNRM);
fig = get(h,'parent');

hd = getappdata(fig,'handles');

imh = getappdata(fig,'imhandle');

if get(hd.ycheck,'Value')
    sl = ceil(get(hd.yslide,'Value')*imdim(2));    
    if ~ishandle(imh(2))
        imh(2) =surf(XDEC(:,sl,:),YDEC(:,sl,:),squeeze(ZDEC(:,sl,:)),squeeze(IMNRM(:,sl,:)));
    else
        set(imh(2),'Xdata',squeeze(XDEC(:,sl,:)),'Ydata',squeeze(YDEC(:,sl,:)),'Zdata',squeeze(ZDEC(:,sl,:)),'Cdata',squeeze(IMNRM(:,sl,:)));
    end
    set(imh(2),'facealpha',get(hd.aslide,'value'));  
  set(hd.ycheck,'String',sprintf('Y: %3i',sl))

elseif ishandle(imh(2)) 
    delete(imh(2));
     set(hd.ycheck,'String',sprintf('Y: %3i',[]))

end
setappdata(fig,'imhandle',imh);
shading flat


%%%%%%%%%%%%%

function redrawz(h,ev,hd)
%Redraw Z slice


global IMNRM ZDEC XDEC YDEC

imdim = size(IMNRM);
fig = get(h,'parent');

hd = getappdata(fig,'handles');

imh = getappdata(fig,'imhandle');

if get(hd.zcheck,'Value')
    sl = ceil(get(hd.zslide,'Value')*imdim(3));    
   
    if ~ishandle(imh(3))
        imh(3) =surf(XDEC(:,:,sl),YDEC(:,:,sl),ZDEC(:,:,sl),IMNRM(:,:,sl));
    else
        set(imh(3),'Xdata',XDEC(:,:,sl),'Ydata',YDEC(:,:,sl),'Zdata',ZDEC(:,:,sl),'Cdata',IMNRM(:,:,sl));
    end

    set(imh(3),'facealpha',get(hd.aslide,'value')); 
    set(hd.zcheck,'String',sprintf('Z: %3i',sl))

elseif ishandle(imh(3)) 
    delete(imh(3));
    set(hd.zcheck,'String',sprintf('Z: %3i',[]))

end
setappdata(fig,'imhandle',imh);
shading flat


%%%%%%%%%%%%%%
function redraw(h,e,hd)

%Redraw all slices

fig = get(h,'parent');
redrawx(h)
redrawy(h)
redrawz(h)

shading flat
colormap gray
% set(imh,'facealpha',.8)
axis(getappdata(fig,'axlim'));

axis equal

%%%%%%
function cleanup(h,e,hd)

%Clear global variables

global IMNRM ZDEC XDEC YDEC

clear global IMNRM ZDEC XDEC YDEC
