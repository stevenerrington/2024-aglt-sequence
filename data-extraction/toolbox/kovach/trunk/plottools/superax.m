function [ax1, ax2, reax, review, refl] = superax(X,Y,varargin)


%[ax1, ax2,reax,review, refl] = superax(ax1,ax2,KEYWORDS) superpimposes 
% ax1 onto ax2 by adjusting axis and view properties so that
%the difference between a set of specified  points is minimized. ax1 and
%ax2 are the axes in the fore and back-ground respectively. reax is the new
%position vector. review is the new view vector and refl has value 1 if the image
% is reflected left-to-right and zero otherwise.   
% 
%   Keywords:
%   
%   'Zoom': when choosing points in X, zoom to box defined by [xmin xmax
%       ymin ymax].
%   'KeepOld': Use marker points from previous call to superax. Specifying Y
%       isn't necessary.
%   'NoRot': no rotation of X.
%   'NoTrans': no translation of X
%   'NoDil':  no dilation of X.
%   'KeepAll': uses all marker points and transformation parameters from
%           previous call to superp.
%   'NoResize': 'Returns error if X and Y are different sizes.    
%   'Transapency': set transparency of ax1.
%   'Colormap': set colormap. Background is converted to true-color image,
%               so isn't affected by colormap setting.
%   'ForegroundPts': Points in the foreground 
%   'BackgroundPts': Points in the background corresponding to points in
%                   the foreground.
%   'PointLabels':
%   'Shading': 'interp' -> bilinear interpolation (default).  'flat' ->
%       nearest point.
%   'Size': followed by factor by which figure heigth and width are reduced
%       from full screen; i.e. .5 makes the figure edges half as big as if it
%       were maximized.

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/plottools/superax.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

%C. Kovach 2002

persistent old_points old_Ly old_fig old_Lx old_ax old_view old_alpha shade old_refl% Ir Jr Id Jd


pos = [0    0.0250    1.0000    0.9158]; %initial position of the figure 

size_reduction_factor = .7;  %reduces figure width by factor 

varg = varargin;

%%%%% Default settings.
zoombox = [];
i = 1;
fig = [];
keep_old_points = 0;
%i = 3;
rot = 1;
dil = 1;
trans = 1;
nl = 0;
keepall = 0;
noresize = 0;
alpha = .6;
cmap = 'jet';
points = {};
Lx = [];
Ly = [];
refl = 0;  %default is no reflection
%%%%%%%%

if ~strcmp(shade,'flat')|~strcmp(shade,'faceted')
    shade = 'interp';
end


  while i <= length(varg)
    
    switch lower(varg{i})
        case 'zoom'
            zoombox = varg{i+1};
        case 'keepold'
            keep_old_points = 1;
            i = i-1;
        case 'norot'
            rot = 0;
            i = i-1;
        case 'nodil'
            dil = 0;
            i = i-1;    
        case 'notrans'
            trans = 0;
            i = i-1;
        case 'keepall'
            keepall = 1;
            keep_old_points = 1;
            i = i-1;
        case 'noresize'
            noresize = 1;
        case 'transparency'
            alpha = varg{i+1};
            old_alpha = alpha;
        case 'colormap'
            cmap = varg{i+1};
        case 'shading'
            shade = lower(varg{i+1});
       case 'foregroundpts'
           Lx = varg{i+1};
       case 'backgroundpts'
           Ly = varg{i+1};
           if isempty(points)
               pts = {1:length(Ly)};
           end
       case 'pointlabels'
           points = varg{i+1};
       case 'figure'
           fig = varg{i+1};
           if ischar(fig)
               fig = figure;
               i = i-1;
           end
       case 'size'
           size_reduction_factor = varg{i+1};
       otherwise
            error(['''',varg{i}, ''' is an invalid keyword'])
    end
    i = i+2;
end
 
 if keepall & ~isempty(old_refl)
    refl = old_refl;
end


if ~ischar(Y) & max(size(Y)) > 1
    
    
    if isempty(fig)
        fig = figure;
    end
    
    if isempty(old_fig)
    
      old_fig = fig;
    elseif isempty(fig)
      fig = figure(old_fig);
    end
    
    ax2 = axes; colormap gray
    s2 = imagesc(Y); axis image                                    %AXis Image
    
    if size(Y,3) < 3
    
        Y = cscale(s2,'rgb');
    end
elseif ~ischar(Y)
    ax2 = Y;
    axes(ax2)
    s2 = get(ax2,'children');
    if isempty(fig)
        fig = gcf;
    end
    old_fig = fig;
    Y = get(s2(end),'cdata');
    if size(Y,3) == 1
        Y = cscale(s2,'rgb');
    end
end


if ~ischar(X) & max(size(X)) > 1
    

    
    figure(fig)
    ax1 = axes;
    i = i+1;
elseif ~ischar(X)
    ax1 = X;
    
    s1 = get(ax1,'children');
    X = get(s1,'CData');
    if strcmp(get(s1,'type'),'surface')
        X = flipud(X);
    end
    if get(ax1,'parent') == fig
        axes(ax1)
    else
 
       
        figure(fig)
        ax1 = axes;
        s1 = pcolor(flipud(X));
        axis image                                                     %Axis  Image
        shading interp, colormap jet
    end
   
end


axes(ax2); cla
image(Y); axis image                                                   %Axis image

        
if gcf~=fig
    figure(fig)
end

X = flipud(X);



set(fig,'units','normalized');

set(fig,'position', pos,'resize','off')

if ~strcmp(shade,'interp')
    X = [X,zeros(size(X,1),1)];
    X = [X;zeros(1,size(X,2))];
    axes(ax1)
    if gcf ~= fig
        figure(fig)
        ax1 = axes
    end
    s1 = pcolor(X);
    shading(shade)
end

i = 1;

pname = 1;

axes(ax2)
axis off

if keep_old_points & ~isempty(old_points)
     Ly = old_Ly;
     points = old_points; 
     ax2_position = axis(ax2);
 elseif isempty(Ly)
  axes(ax2);
  
  if gcf ~= fig
      figure(fig)
      ax2 = axes;
      s2 = image(Y);
  end
  vw = get(ax1,'view');
  if vw(1) == 0 
  end
  
  while ~isempty(pname)
    [gui_out,gui_handle] = getlabel;
    delete(gui_handle)
    zoom off
    %pname = input(['\n>Label for point ',num2str(i),' (<Ret> = end): '],'s');
    pname = gui_out;
    if isempty(pname), break, end
    figure(fig)
    Ly = [Ly;getpoint(ax2)];    
       
    points = cat(2,points,{pname});
    
    i = i+1;
  end
  ax2_position = axis(ax2);
   
end

old_points = points;
old_Ly = Ly;


if keepall & ~isempty(old_Lx)
    Lx = old_Lx;
    axes(ax1) 
    s1 = pcolor(X); shading(shade); 
    if refl
        X = reflect(X,Lx,s1);
    end
    colormap(cmap), axis image                             %axis image
elseif isempty(Lx)
    vw = get(ax1,'view');
    axes(ax1),  set(ax1,'view',[vw(1) 90])
    figure(fig)
    s1 = pcolor(X);
    shading(shade), colormap(cmap)
    axis image % colormap gray 

  if ~isempty(zoombox)
    axis(zoombox)
  end

  for j = points
    
    %fprintf(1,['\n>Click on point corresponding to ',char(j)])
    title_handle = title(['Click on point corresponding to ',char(j)],'fontsize',15);
    Lx = [Lx;getpoint(ax1)];
  end
  delete(title_handle)

  refl = isflipped(Lx,Ly);
  
  if refl
    [X,Lx] = reflect(X,Lx,s1);
  end
  
  old_refl = refl;
  
end






axis manual
old_Lx = Lx;

size_reduction = [1 1 size_reduction_factor*[1 1]];
set(fig,'position',pos.*size_reduction)            % Sets size 

mx = mean(Lx); %centroid of L1
m = mean(Ly); %centroid of L2



if ~keepall | isempty(old_ax) | isempty(old_view) %~isempty(Ir)&~~isempty(Jr)&~isempty(Id)&~isempty(Jd)
   
   reax = get(ax1,'position');
   review = get(ax1,'view');
   
  if dil
    [reax,L] = dial(ax1,Lx,Ly);                           %Local function, Dilation of Xm.  
  end
 
   if trans
    [reax,L] = translate(ax1,L,Ly);                      %Local function translation to give same mean L.
    % reax(1:2) = axtrs; 
   end


  if rot
    [review,L] = rotate(ax1,L,Ly);                     %Local function, Rotation of Xm.
  end
  
  old_ax = get(ax1,'position');
  old_view = get(ax1,'view');
else     
  axes(ax1)
  reax = old_ax;
  review = old_view;
  set(ax1,'position',reax);
  set(ax1,'view',review);
  shading(shade)
  colormap(cmap)
end

if keep_old_points & ~isempty(old_alpha)
    alpha = old_alpha;      %Sets transparency to old transparency
end

%close(gcf)

axes(ax2)
axis(ax2_position)
axes(ax1), axis off
%set(fig,'renderer','painter')
set(s1,'facealpha',alpha)
reax = get(ax1,'position');
review = get(ax1,'view');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reax,L] = translate(ax,Lx,Ly);    %finds appropriate position of lower left corner

mx = mean(Lx); %centroid of L1
m = mean(Ly); %centroid of L2

dm = mx-m;

reax = get(ax,'position');

L = Lx - dm(ones(length(Lx),1),:);

reax(1:2) = reax(1:2) -dm;              %(end:-1:1);  %Flipped dm 
set(ax,'position',reax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [review,L] = rotate(ax,Lx,Ly)  %finds appropriate 'view' property


mx = mean(Lx); %centroid of Lx
m = mean(Ly); %centroid of Ly

Lxc = Lx - mx(ones(size(Lx,1),1),:);
Lyc = Ly - m(ones(size(Ly,1),1),:);

fig = get(ax,'parent');
fig_units = get(fig,'units');
set(fig,'units','pixels');
fig_pos = get(gcf,'position');

Lxc(:,1) = Lxc(:,1)*fig_pos(3);     %Lx and Ly are in units normalized to figure height and width,
Lxc(:,2) = Lxc(:,2)*fig_pos(4);     %so need to scale according to pixel height and width of figure.
Lyc(:,1) = Lyc(:,1)*fig_pos(3);
Lyc(:,2) = Lyc(:,2)*fig_pos(4);

inv_norms = (sum(Lxc.^2,2).*sum(Lyc.^2,2)).^-.5;
sin_thetas = (Lxc(:,1).*Lyc(:,2) - Lyc(:,1).*Lxc(:,2)).*inv_norms;  %cross product
cos_thetas = sum(Lxc.*Lyc,2).*inv_norms;                            %inner product

TH = real(j*log(cos_thetas - j*sin_thetas));               %angle of rotation

%TH = asin((Lxc(:,1).*Lyc(:,2) - Lyc(:,1).*Lxc(:,2)).*(sum(Lxc.^2,2).*sum(Lyc.^2,2)).^-.5); 

wgt = sqrt(sum((Lx-mx(ones(size(Lx,1),1),:)).^2,2));
wgt = wgt./sum(wgt);

%length_Lxc = sqrt(sum(Lxc.^2,2));
%length_Lyc = sqrt(sum(Lyc.^2,2));

vw = get(ax,'view');
axs = get(ax,'position');


thrad = -sign(vw(2))*(TH'*wgt);                              %average rotation, weighted by distance from center of rotation (to reduce the effect of error).
th = thrad*180/pi;                                          %negative if elevation of view is 90.
%adj = [ax(4)./(ax(4)*abs(cos(th))+ax(3)*abs(sin(th))), ax(3)./(ax(4)*abs(sin(th))+ax(3)*abs(cos(th)))];
Lc = ([cos(thrad) -sin(thrad); sin(thrad) cos(thrad)]*Lxc')';  %.*adj(ones(1,length(Lx)),:);
L = Lc + m(ones(size(Lc,1),1),:);

review = vw;
review(1) = review(1) + th;

reax = axs;
reax(3:4) = [abs(axs(3)*cos(thrad)) + abs(axs(4)*sin(thrad)), abs(axs(3)*sin(thrad)) + abs(axs(4)*cos(thrad))];
reax(1:2) = reax(1:2) + (axs(3:4)-reax(3:4))./2;
set(ax,'position',reax)
set(ax,'view',review)


set(fig,'units',fig_units)

L = L./fig_pos(ones(length(L),1),3:4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reax,L] = dial(ax,Lx,Ly)      %finds appropriate scaling for axis width and heighth


fig = get(ax,'parent');
fig_units = get(fig,'units');
set(fig,'units','pixels');
fig_pos = get(gcf,'position');

Lx(:,1) = Lx(:,1)*fig_pos(3);     %Lx and Ly are in units normalized to figure height and width,
Lx(:,2) = Lx(:,2)*fig_pos(4);     %so need to scale according to pixel height and width of figure.
Ly(:,1) = Ly(:,1)*fig_pos(3);
Ly(:,2) = Ly(:,2)*fig_pos(4);

mx = mean(Lx); %centroid of Lx
m = mean(Ly); %centroid of Ly

Dy = sqrt(sum((Ly - m(ones(size(Lx,1),1),:)).^2,2));
Dx = sqrt(sum((Lx - mx(ones(size(Lx,1),1),:)).^2,2));
wgt = Dy./sum(Dy);
 

s = sum((Dy./Dx).*wgt);
s = [s,s];


axs = get(ax,'position');
reax = axs;
reax(3:4) = axs(3:4).*s;                            %(end:-1:1);    flipped
%reax(1:2) = reax(1:2) + (axs(3:4)-reax(3:4))./2;

Lx = [Lx;reax(1:2).*fig_pos(3:4)];                                     %(2:-1:1)];         "
L = (Lx - mx(ones(size(Lx,1),1),:)).*s(ones(size(Lx,1),1),:) + mx(ones(size(Lx,1),1),:);
L = L./fig_pos(ones(length(L),1),3:4);
reax = [L(end,:),reax(3:4)];                        %flipped: colon replaces '2:-1:1'
L = L(1:end-1,:);
set(ax,'position',reax)
set(fig,'units',fig_units)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = getpoint(ax)
axes(ax)
k = 1;
ptr = get(gcf,'pointer');
units = get(gcf,'units');
set(gcf,'units','normalized')
set(gcf,'pointer','fullcrosshair')

while k
  axes(gca)  
  k = waitforbuttonpress;
  P = get(gcf,'currentpoint');
end

set(gcf,'pointer',ptr)
set(gcf,'units',units) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  f = isflipped(Lx,Ly)

%Detects wether Lx is flipped wrt to Ly

mx = mean(Lx);
my  =mean(Ly);
Dfx = diff(Lx);                                                                         %
Dfy = diff(Ly);                                                                         %
Xmidpts = [Lx(1:end-1,:) + .5*Dfx - mx(ones(length(Dfx),1),:), zeros(length(Dfx),1)];   %   Using cross-product to find if 
Ymidpts = [Ly(1:end-1,:) + .5*Dfy - my(ones(length(Dfy),1),:), zeros(length(Dfy),1)];   %   plot is reflected with respect to
Dfx = [Dfx,zeros(length(Dfx),1)];                                                       %   background (i.e. flipped left-to-
Dfy = [Dfy,zeros(length(Dfy),1)];                                                       %   right, or up to down).
t = sum(cross(Dfx',Xmidpts').*cross(Dfy',Ymidpts'),2);                                  %
f = t(3) < 0;

%if t(3) < 0                                               %  
%    vw = get(ax,'view');                                  %
%    vw(2) = vw(2)-180;                                    %  possible bug when th = 0;
%    set(ax,'view',vw);                                    %
%end                                                       %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RX,RLx]  = reflect(X, Lx, s)

%Returns X reflected about midpoint of figure.

RX = fliplr(X);

ax = get(s,'parent');
pos = get(ax,'position');
shade = get(s,'FaceColor');

midpt = [pos(1) + .5*pos(3), pos(2) + .5*pos(4)];

if ~strcmp(shade,'interp')
    RX = [RX(:,2:end),X(:,end)];
end
    
RLx = Lx;

RLx(:,1) = 2*midpt(ones(length(Lx),1),1) - Lx(:,1);    %reflection
 

set(s,'CData', RX)