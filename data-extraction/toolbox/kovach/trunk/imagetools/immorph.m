function Mstruc = immorph(X,varargin)

% 
% IMMORPH morphs images by applying a piecewise linear transformation over a delaunay tesselation.
%   The image interpolation is carried out by GRIDDATA
%
% How to use:
% 
% 1. load two images (X and Y) into Matlab  using imread, eg. X = imread('imageName.ext','ext').
% 2. Convert images to double and scale to [0 1], eg  X = double(X)/255.
% 3. Add fiducial points with the immorph command, M = immorph(X,Y)
%         
%         A. Right clicking will add a new active point to each image. The active point is indicated with a cross.
%         B. Left clicking will move an active point to a new position.
%         C. Right clicking again will de-activate the current active point.
%         D. Left clicking over an inactive point will activate it if no other point is already active.
% 
% 4. Press the space bar when finished adding fiducial points.
% 5. IMMOPRH returns a data structure, M, with the following fields
%         
%             X: matrix with the original image X
%             Y: matrix with the original image Y
%             VX: fiducial point locations on X
%             VY: fiducial point locations on Y
%             Tess: matrix containing indices of the ponts in VX and VY that form the vertices of each tesselation.
%                 Tesselation is done with the DELAUNAY function.
%             mapto: mapping of each pixel in X to a pixel in Y
%             makepic: function handle to the function that generates the morph.
%             tessmap: function handle to a function that generates a map of tesselations.
% 
% 6. IMMORPH takes the following optional keywords
% 
%     'align': Aligns Y with X by rotating and scaling Y to minimize squared difference between the fiducial points.
%     'vertices': Followed by VX and VY. Adds vertex (fiducial) points from the matrices VX and VY:  
%         M = immorph(X,Y,'vertices',VX,VY)  adds points in VX to X and corresponding points VY to Y.These can be edited subsequently.
%     'noget': Runs IMMORPH without acquiring additional fiducial points. Uses only points supplied with the
%         'vertices' option.
%     'fix': Followed by 1 or 2. Disallows any changes to points in X ( if 1 ) or Y ( if 2 ). No additional points can be
%         added to either image if either of them is fixed.
%         
%     
% Generating Morphs:
% 
%    1. Morphs are generated with the M.makepic function handle, which is called using FEVAL:
%      The structure M, itself, has to be passed as an argument to M.makepic. 
% 
%         morphPic = feval(M.makepic,M,morphlevel).
%      
%       where morphlevel is the fractional change between X and Y (morphlevel = .5 gives a 50%
%           morph).
%
%     2. Morphing level can be specified independently for each pixel by using the keyword 'mask' followed by 
%        a matrix the same size as the image: 
%               
%        morphPic = feval(M.makepic,M,morphlevel,'mask',maskMatrix)
% 
%     3. Shading and warping can be controlled independently using the keywords 'shadingscale' and 'warpingscale' respectively.
%         In both cases 0 preserves X and 1 preserves Y. To warp X to Y use:
%         
%             morphPic = feval(M.makepic,M,0,'warpingscale',1)  or 
%             morphPic = feval(M.makepic,M,1,'shadingscale',0)

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/imagetools/immorph.m $
% $Revision: 405 $
% $Date: 2013-11-20 18:22:04 -0600 (Wed, 20 Nov 2013) $
% $Author: ckovach $
% ------------------------------------------------

% Written by Christopher Kovach 2006,2012
            
       
VX = [];
VY = [];
if isstruct(X)
    VX = X.VX;
    VY = X.VY;
    Y = X.Y;
    X = X.X;
else
    Y = varargin{1};
    varargin = varargin(2:end);
end


% eraseUnmappedPixels = 1;
getpoints = 1;
roundmapto = 1;
fixax = 0; %excludes vertices in one axis from editing. 1 for left, 2 for right.
align = 0;
% Tesswgt = .5; %Tesselation is drawn on a map that is interpolated between VX (0) and VY (1) according to this weight.
Tesswgt = .5; 
i = 1;
while i <= length(varargin)
    switch lower(varargin{i})
        case 'vertices'
            VX = varargin{i+1};
            VY = varargin{i+2};
            i = i+2;
         case 'noget',
            getpoints  = 0;
%          case 'noerase',
%             eraseUnmappedPixels = 0;
         case 'noround'
              roundmapto = 0;
         case 'fix'
               fixax = varargin{i + 1};
               i = i+1;
        case 'align'
            align = varargin{i + 1};
               i = i+1;
        case 'tess'
           Tesswgt = varargin{i + 1};
               i = i+1;
        otherwise
            error('%s is not a valid keyword',varargin{i})
    end
    
    i = i+1;
end

tempsize = [max([size(X,1),size(Y,1)]),max([size(X,2),size(Y,2)]),max([size(X,3),size(Y,3)])];
template = zeros(tempsize);
if ~isequal(size(X),tempsize)
   T = template;
   T(1:size(X,1),1:size(X,2),:) = X;
   X = T;
end

if ~isequal(size(Y),tempsize)
   T = template;
   T(1:size(Y,1),1:size(Y,2),:) = Y;
   Y = T;
end

clear T



VXin = VX;
VYin = VY;

strs = {' has','s have'};
str = strs{(abs(size(VXin,1)-size(VYin,1)) > 1)+1};

if size(VXin,1)>size(VYin,1)
    
    warning('Number of Vertex points doesn''t match. %i point%s been added to 2nd image',size(VXin,1)-size(VYin,1),str)
    VYin(end+1:size(VXin,1),:) = VXin(size(VYin,1)+1:end,:);

elseif size(VXin,1)<size(VYin,1)
    
    warning('Number of Vertex points doesn''t match. %i point%s been added to 1st image',size(VYin,1)-size(VXin,1),str) %#ok<*WNTAG>
    VXin(end+1:size(VYin,1),:) = VYin(size(VXin,1)+1:end,:);
    
end    

szeX = size(X);
% szeY = size(Y);


if getpoints
    [VX,VY] = getvertices(X,Y,VXin,VYin,fixax,Tesswgt);
else
    VX = VXin;
    VY = VYin;
end

if ~isscalar(align)
    [Y,VY] = alignImages(Y,align,VY);    
elseif align
    [Y,VY] = alignImages(Y,VX,VY);
end

fprintf('Computing pixel map... |')

VTess = (VX*(1 - Tesswgt) + VY*Tesswgt);
Tess = delaunay(VTess(:,1),VTess(:,2)); % delaunay tesselation of fiducial points;

[pixmapI,pixmapJ] = ind2sub(szeX(1:2),1:prod(szeX(1:2)));
pixmap = [pixmapI;pixmapJ];
% [pixmapYI,pixmapYJ] = ind2sub(szeY(1:2),1:prod(szeY(1:2)));

mapto = zeros(2,size(pixmapI,2));
%insidemapY = zeros(szeY(1:2));
insidemapX = zeros(szeX(1:2));
rotor = '/-\|';
detA = 0;
for i = 1:size(Tess,1)
    
    fprintf('\b%s',rotor(mod(i,4)+1))
    
    %find points within each tesselation triangle
    insidepts = find(inpolygon(pixmapJ,pixmapI,VX(Tess(i,:),1),VX(Tess(i,:),2)));
    %insideptsY = inpolygon(pixmapYI,pixmapYJ,VY(Tess(i,:),1),VY(Tess(i,:),2));
    
%    insidemapY(pixmapYI(find(insideptsY)),pixmapYJ(find(insideptsY))) = i;
%    insidemapX(pixmapI(find(insideptsX)),pixmapJ(find(insideptsX))) = i;

    %insidemapY(find(insideptsY)) = i;
    insidemapX(insidepts) = i;
 

    %Linear transform: y = Ax + b; Since vy = A*(vx - vx0) + vy0, and [vy1 vy2] = A*([vx1,vx2] - [vx0,vx0]) + [vy0,vy0]
    % then A = ([vy1 vy2] - [vy0,vy0]) * ([vx1,vx2] - [vx0,vx0])^-1;
    A = ([VY(Tess(i,2),:)',VY(Tess(i,3),:)'] - [VY(Tess(i,1),:)',VY(Tess(i,1),:)']) * ([VX(Tess(i,2),:)',VX(Tess(i,3),:)'] - [VX(Tess(i,1),:)',VX(Tess(i,1),:)'])^-1;
   
    %Vertices are in XY form, hence flipud to make mapto IJ.
    mapto(:,insidepts) =  flipud(A * ([pixmapJ(insidepts);pixmapI(insidepts)] - VX(Tess(i,1),:)'*ones(1,length(insidepts))) + VY(Tess(i,1),:)'*ones(1,length(insidepts))) ;
    detA(i)=det(A);
end

%if ~eraseUnmappedPixels
 if roundmapto
      mapto = round(mapto);
  end
    

  mapto(mapto <= 0) = pixmap(mapto <= 0);
  mapto(1,(mapto(1,:) > size(Y,1))) = pixmap((mapto(1,:) > size(Y,1)));
  mapto(2,(mapto(2,:) > size(Y,2))) = pixmap((mapto(2,:) > size(Y,2)));
  
  [Iind,Jind] = ind2sub(size(Y),find(isnan(mapto(1,:))));
  mapto(:,isnan(mapto(1,:))) = [Iind; Jind];
  
  %else
  %unmappedX = pixmap(mapto == 0);  
  %mappedY = mapto(mapto(
  %X(unmapped(1,:),unmapped(2,:)) = 0;
  %Y(ismember(

Mstruc.X = X;
Mstruc.Y = Y;
Mstruc.VX = VX;
Mstruc.VY = VY;
Mstruc.Tess = Tess;
Mstruc.detA=detA;
Mstruc.detA=Tesswgt;
if roundmapto
  Mstruc.mapto = sub2ind(size(X),mapto(1,:),mapto(2,:));
else
  Mstruc.mapto = mapto;
end

Mstruc.makepic = @makepicfun;
Mstruc.tessmap = @tessmapfun;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M,MedVs] = makepicfun( Mstruc , pcnts , varargin ) 

% Function to generate morphed image. Mstruc is the morph image structure
% returned by immorph.

shadingscale = .5; % 0 Preserves shading from 1st picture and 1 preserves from 2nd.
warpingscale = .5; % "      "    warping    "       "       "       "       "
mask = .5;          %mask can be a matrix identifying the morphing degree for each pixel.

block_neg_det = true; %%%  regard a negative determinant as “behind”
foreground = 2;    %%%  in the event of a negative determinant, which image is in front?

% method = 'linear';
i = 1;
OthogonalWarping = 0;
gray = 0;
align = false;
while i <= length(varargin)
    switch lower(varargin{i})
        case 'shadingscale'
            shadingscale = varargin{i + 1};
            i = i+1;
%         case 'method'
%             method = varargin{i + 1};
%             i = i+1;
         case 'warpingscale'
            warpingscale = varargin{i + 1};
            i = i+1;
        case 'orthwarp'   %Warps images perpendicular to the difference. 
            OthogonalWarping = 1;
            orthscale = varargin{i + 1};
            i = i+1;

        case 'gray'
            gray = 1;
        case 'foreground'
            foreground =  varargin{i + 1};
            i = i+1;
        case 'mask'
            mask = varargin{i + 1};
            i = i+1;
        case 'align'
            align =  varargin{i + 1};
            i = i+1;
        otherwise
            error('%s is not a valid keyword',varargin{i})
    end
    i = i+1;
end



XX = Mstruc.X;
YY = Mstruc.Y;
VX = Mstruc.VX;
VY = Mstruc.VY;

if ~isscalar(align)
    [YY,VY] = alignImages(YY,align,VY);    
elseif align
    [YY,VY] = alignImages(YY,VX,VY);
end 

if size(Mstruc.mapto,1) == 1
  [I,J] = ind2sub(size(XX),Mstruc.mapto);
  mapto = [I;J];
else
  mapto = Mstruc.mapto;
end

if gray 
    XX = mean(XX,3);
    YY = mean(YY,3);
end
    
szeX = size(XX);
szeY = size(YY);

[pixmapI,pixmapJ] = ind2sub(szeX(1:2),1:prod(szeX(1:2)));
% pixmap = [pixmapI;pixmapJ];


if any(Mstruc.detA < 0) && block_neg_det
   [xx,tind] = Mstruc.tessmap(Mstruc,pcnts + randn*.00001);  %%% Joggle slightly 
%    dA = 1-pcnts + pcnts*Mstruc.detA ;
%    detval = tind*diag(dA);
   
   nd = any(sum(tind,2)>1,2);
   if any(nd)
       if length(shadingscale)==1
           shadingscale=shadingscale*ones(szeX);
       end
       shadingscale(nd) = foreground-1;
       
         shadingscale = ones(2,1)* shadingscale(:)';

   end
    
end

if length(mask) ~= 1
    
    mask = ones(2,1)*mask(:)';

end    
M = {};
MedVs = {};
j = 0;
fprintf('\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t')
for s = pcnts
   
   j = j+1;
   fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') 
   fprintf('%2.0d of %2.0d ( %3.0f%% morphing )',j,length(pcnts),s*100)
    
   MedVs{end+1} = VX + s*(VY - VX);
   for cl = 1:size(XX,3)
    X = XX(:,:,cl);
    Y = YY(:,:,cl);
    sizem = ceil(szeX(1:2) + s * (szeY(1:2) - szeX(1:2)));
    
%     [pixmapMI,pixmapMJ] = ind2sub(sizem,1:prod(sizem));
    
    warps = s;
    shades = s;
    
    if warpingscale == 0
        warps = 0;
    elseif warpingscale == 1
        warps = 1;
    end
    
    if shadingscale == 0
        shades = 0;
    elseif shadingscale == 1
        shades = 1;
    end
    
    
      MedPs = ((1 - mask )*( 1 - warps ) * ( 1 - warpingscale ) .* [pixmapI;pixmapJ] + mask * warps * warpingscale .* mapto)./( mask .* warps * warpingscale + (1 - mask).*(1 - warps) * (1 - warpingscale) );
      MedZ =  ((1 - mask(1,:)' )*( 1 - shades ) * ( 1 - shadingscale(1,:)' ) .* X(:) + mask(1,:)' * shadingscale(1,:)' * shades .* Y(sub2ind(szeY(1:2),round(mapto(1,:)),round(mapto(2,:)))'))./(mask(1,:)' .* shades * shadingscale(1,:)' + (1 - mask(1,:)').*(1 - shades) * (1 - shadingscale(1,:)') );
      
    if OthogonalWarping    %Warping orthogonal to difference vector of pixel positions
        orthvec = [0 -1 ; 1 0]*(mapto-[pixmapI;pixmapJ]) *orthscale;
        MedPs = MedPs + orthvec;    
    end
     
    [mshx,mshy] = meshgrid(1:sizem(2),1:sizem(1)); 
    
    warning off MATLAB:griddata:DuplicateDataPoints
     
%     MdX = reshape(MedPs(2,:),szeX(1:2));
%     MdY = reshape(MedPs(1,:),szeX(1:2));
%     MdZ = reshape(MedZ,szeX(1:2));
%     
%     mz(:,:) = griddata(MdX,MdY,MdZ,mshx,mshy,method,{'QJ'});  
%%% New fcn as of 2011 is more efficient
    interpfun = TriScatteredInterp(MedPs(2,:)',MedPs(1,:)',MedZ);
    mz = reshape(interpfun(mshx(:),mshy(:)),szeX(1:2));
    
    badpix = find(isnan(mz));
    
    if ~isempty(badpix)
        
%         [bdi,bdj] = ind2sub(sizem,badpix);
%         
%         for bdind = 1:length(bdi) 
%             bdis = max(bdi(bdind)-1),1):min(bdi(bdind)+1,sizem(1)); 
%             bdjs = max(bdj(bdind)-1),1):min(bdj(bdind)+1,sizem(1)); 
%             bdpix = meshgrid(
%          try
%             pixmeans = diag(mz(bdi+1,bdj) + mz(bdi-1,bdj) + mz(bdi,bdj-1) + mz(bdi,bdj+1) + mz(bdi-1,bdj-1) + mz(bdi+1,bdj+1)...
%                  + mz(bdi+1,bdj-1) + mz(bdi-1,bdj+1))/8;
%             mz(sub2ind(sizem,badpix)) = pixmeans;
%         end        
        
        %badpix = find(isnan(mz));
        mz(badpix) = (( 1 - shades ) * ( 1 - shadingscale(1) ) * X(badpix) + shadingscale(1) * shades * Y(badpix))./(shades * shadingscale (1)+ (1 - shades) * (1 - shadingscale(1)) );
        
    end  
   
    
    MZ(:,:,cl) = mz;
    
    % MZ(:,cl) = fastinterp(MdX,MdY
    %for j = 1:size(Tess,1)
    
      %  insideptsM = find(inpolygon(pixmapMJ,pixmapMI,MedVs(Tess(j,:),1),MedVs(Tess(j,:),2)));
        %interpolate points on the 
%        mz = interp2(MedPs(1,insidemapX == j),MedPs(2,insidemapX == j),X(insidemapX == j) + s * (X(insidemapX == j) - Y(sub2ind(szeY(1:2),round(mapto(1,insidemapX == j)),round(mapto(2,insidemapX == j))))') , pixmapMI(insideptsM),pixmapMJ(insideptsM)); 
        %m(pixmapMI(insideptsM),pixmapMJ(insideptsM)) = mz;
       
        
        %end
    end   
     M{end+1} = MZ;
end

if length(M) == 1
    M = M{1};
    MedVs = MedVs{1};
end 









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VX,VY] = getvertices(X,Y,VX,VY,fixax,tesswgt)

PtSelectionTol = 10;

fig = figure; 
colormap gray
ax1 = subplot(1,2,1); imagesc(X);
axis image, axis off, hold on
ax2 = subplot(1,2,2); imagesc(Y);
axis image, axis off, hold on
% zoom on
% z = zoom(fig);
% set(z,'ButtonDownFilter',@(a,b)false)
% set(z,'ActionPreCallback',@(a,b)set(a,'userdata',0))
setappdata(fig,'tesswgt',tesswgt)
axvec = [ax1 ax2]; 

colors = {'b','g','r','c','m','y','w'};

%VX = [];
%VY = [];
p1 = [];
p2 = [];
j = 0;

quit = 0;
plotcolr = [];

    plotopts= {'markersize',12};

    pass2plot1 = {};
   pass2plot2 = {};
    
   for k = 1:size(VX,1)
%         
%         if k == currvertind
%             add = '+';
%         else
            add = '.';
%         end
%         
        pass2plot1{1} =  VX(k,1);
        pass2plot1{2} =  VX(k,2);
        pass2plot1{3} =  [colors{mod(k,length(colors))+1},add];
        
        pass2plot2{1} =  VY(k,1);
        pass2plot2{2} =  VY(k,2);
        pass2plot2{3} =  [colors{mod(k,length(colors))+1},add];

        p1(end+1) = plot(ax1,pass2plot1{:},plotopts{:});
        p2(end+1) = plot(ax2,pass2plot2{:},plotopts{:});

   end
   
   
u = uicontrol('visible','on','units','normalized','position',[.46 .05 .1 .05],...
   'style','togglebutton','string','tesselate');
setappdata(fig,'buttonh',u);
setappdata(fig,'axes',[ax1 ax2]);
set(u,'Callback',@trplot)
setappdata(fig,'VX',VX)
 setappdata(fig,'VY',VY)
    
while ~quit
    
    colorind = mod(j,length(colors))+1;
    
    figure(fig)    
    refresh(fig)
    pause(.05)
    waitforbuttonpress
%     waitfor(fig,'userdata',0)
%     drawnow
    %     pause(1)
   
            
    %get(gcf,'selectiontype')
    if isempty(get(gcf,'currentcharacter'))
        quit = 0;
    elseif isspace(get(gcf,'currentcharacter'))
        break
    end
    
    if gca == ax1
%         currax = ax1;
        currV = VX;
    elseif gca == ax2
%         currax = ax2;
        currV = VY;
    else
%         currax = [];
%         currV = [];
        continue
    end

    currpt = get(gca,'currentpoint'); 
    currpt = currpt(1,1:2);
    
    if strcmp(get(gcf,'selectiontype'),'normal')
            
      try
        dv = ones(size(currV,1),1)*currpt - currV;
        dv = sqrt(sum(dv.^2,2));
        [mn,currvertind] = min(dv);
        if mn < PtSelectionTol
          if gca == ax1
              pa = p1;
              pb = p2;
          elseif gca == ax2
              pa = p2;
              pb = p1;
          end
          set(pa(currvertind),'xdata',currV(currvertind,1),'ydata',currV(currvertind,2),'marker','+')
          set(pb(currvertind),'marker','+')
          refresh(fig)
        else
          continue
        end    
      catch %#ok<*CTCH>
          continue
      end
    elseif strcmp(get(gcf,'selectiontype'),'alt')
       
       if fixax ~= 0, continue, end
       
 
       if fixax ~= 0, continue,end
 
       try 
        curraxdim = axis(gca);
        
        if ~((currpt(1) < 0) || (currpt(2) < 0) || (currpt(1) > curraxdim(2)) || (currpt(2) > curraxdim(4)))
         currV = [currV;round(currpt(1,:))]; %#ok<*AGROW>
         currvertind = size(currV,1);
         VX = [VX;round(currpt(1,:))];
         VY = [VY;round(currpt(1,:))];
         plotcolr(end+1) = colorind;
%          axes(ax1)
         p1(end+1) = plot(ax1,VX(end,1),VX(end,2),[colors{colorind},'+'],plotopts{:});
%          axes(ax2)
         p2(end+1) = plot(ax2,VY(end,1),VY(end,2),[colors{colorind},'+'],plotopts{:});
         refresh(fig)
         j = j+1;
        else
         continue
        end
       catch
           continue
       end
    elseif strcmp(get(gcf,'selectiontype'),'extend')
       
       if fixax ~= 0, continue,end
       
       try
        dv = ones(size(currV,1),1)*currpt - currV;
        dv = sqrt(sum(dv.^2,2));
        [mn,currvertind] = min(dv);
       if mn < PtSelectionTol
          
%         currV(currvertind,:) = [];
        VX(currvertind,:) = [];
        VY(currvertind,:) = [];
       
        delete(p1(currvertind))
        delete(p2(currvertind))
        
        p1(currvertind) = [];
        p2(currvertind) = [];
        plotcolr(currvertind) = [];
        
        currvertind = [];
        
        refresh(fig) 
       else
           continue
       end
       catch caterr %#ok<NASGU>
           continue
       end
    else
%         currvertind = [];
	 	continue;   
    end
    
    if ~isempty(p1), set(p1,'zdata',1),end
    if ~isempty(p2), set(p2,'zdata',1),end
    
   
    setappdata(fig,'VX',VX)
    setappdata(fig,'VY',VY)
    
%     currvert = currV(currvertind,:);
    
    %add = '';
%     for k = 1:size(currV,1)
%         
%         if k == currvertind
%             add = '+';
%         else
%             add = '.';
%         end
%         
%         pass2plot1{3*(k-1)+1} =  VX(k,1);
%         pass2plot1{3*(k-1)+2} =  VX(k,2);
%         pass2plot1{3*(k-1)+3} =  [colors{plotcolr(k)},add];
%         
%         pass2plot2{3*(k-1)+1} =  VY(k,1);
%         pass2plot2{3*(k-1)+2} =  VY(k,2);
%         pass2plot2{3*(k-1)+3} =  [colors{plotcolr(k)},add];
%               
%     end
    
    set(gcf,'selectiontype','normal')
    
%     delete(p1)
%     axes(ax1)
%     p1 = plot(pass2plot1{:});
%         
%     delete(p2)
%     axes(ax2)
%     p2 = plot(pass2plot2{:});

    while 1
        
       isKey = waitforbuttonpress;
       %get(gcf,'selectiontype')
       
       if ~strcmp(get(gcf,'selectiontype'),'normal') || isKey
            break
       end
       
       if fixax ~= 0
         if axvec(fixax) == gca, continue, end
       end
       
       currpt = get(gca,'currentpoint');
       currpt = round(currpt(1,1:2));
       curraxdim = axis(gca);
       
       if ~((currpt(1) < 0) || (currpt(2) < 0) || (currpt(1) > curraxdim(2)) || (currpt(2) > curraxdim(4)))
         if gca == ax1
            VX(currvertind,:) = currpt;
            set(p1(currvertind),'xdata',VX(currvertind,1),'ydata',VX(currvertind,2))
        elseif gca == ax2
            %currpt = get(gca,'currentpoint');
            %currpt = round(currpt(1,1:2));
            VY(currvertind,:) = currpt;
            set(p2(currvertind),'xdata',VY(currvertind,1),'ydata',VY(currvertind,2))
        end
      end
        %refresh(fig)
       
    end
    
     set(p1(currvertind),'marker','.')
     set(p2(currvertind),'marker','.')
     
    setappdata(fig,'VX',VX)
    setappdata(fig,'VY',VY)
    
    trh = getappdata(fig,'trh');
    if get(u,'value') 
        trplot(u)
    elseif all(ishandle(trh(:)))
        set(trh,'visible','off')
    end
end    

close(fig)
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trplot(bh,ev,hd) %#ok<INUSD>

h = get(bh,'parent');
trh = getappdata(h,'trh');
if ~isempty(trh)&&all(ishandle(trh(:))), delete(trh), end
if get(getappdata(h,'buttonh'),'Value')==0,
    return,
end
VX = getappdata(h,'VX');
VY = getappdata(h,'VY');
axs = getappdata(h,'axes');
wgt = getappdata(h,'tesswgt');


if isempty(VX), return, end
VM = (1-wgt)*VX+wgt*VY;
D=DelaunayTri(VM);
trh = [];
axes(axs(1)) %#ok<MAXES>
trh(:,1) = triplot(D.Triangulation,VX(:,1),VX(:,2));
axes(axs(2)) %#ok<MAXES>
trh(:,2) = triplot(D.Triangulation,VY(:,1),VY(:,2));

setappdata(h,'trh',trh)
% setappdata('VX',VX)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TM,Ind] = tessmapfun(Mstruc,pcnt)

% A function that returns a map of pixels showing the index of the
% tesselation to which they belong.

szeX = size(Mstruc.X);
% szeY = size(Mstruc.Y);

Tess = Mstruc.Tess;
VX = Mstruc.VX;
VY = Mstruc.VY;

Ind = logical(spalloc(prod(szeX(1:2)),length(Tess),prod(szeX(1:2))));
%Ind = logical(spalloc(prod(szeX),length(Tess),prod(szeX)));

MDX = VX + pcnt*(VY - VX);
[pixmapXI,pixmapXJ] = ind2sub(szeX(1:2),1:prod(szeX(1:2)));
% pixmapX = [pixmapXI;pixmapXJ];

TM = zeros(szeX(1:2));

% if picnum == 1
%     [pixmapXI,pixmapXJ] = ind2sub(szeX(1:2),1:prod(szeX(1:2)));
%     pixmapX = [pixmapXI;pixmapXJ];
%     
%     %Js = 1:szeX(2);
%     %Is = 1:szeX(1);
%     TM = zeros(szeX(1:2));
% elseif picnum == 2
%     [pixmapYI,pixmapYJ] = ind2sub(szeY(1:2),1:prod(szeY(1:2)));
%     pixmapY = [pixmapYI;pixmapYJ];
%     
%     %Js = 1:szeY(2);
%     %Is = 1:szeY(1);
%     TM = zeros(szeX(1:2));
% else
%     error('Picture number must be 1 or 2')
% end


for i = 1:size(Tess,1)
    
    
    %find points within each tesselation triangle
%   if picnum == 1
%     insideptsX = find(inpolygon(pixmapXJ,pixmapXI,VX(Tess(i,:),1),VX(Tess(i,:),2)));
%     TM(insideptsX) = i;    
%   elseif picnum == 2
%     insideptsY = find(inpolygon(pixmapYJ,pixmapYI,VY(Tess(i,:),1),VY(Tess(i,:),2)));
%     TM(insideptsY) = i;
%   end

    insideptsX = inpolygon(pixmapXJ,pixmapXI,MDX(Tess(i,:),1),MDX(Tess(i,:),2));
    TM(insideptsX) = i;    
    Ind(:,i) = insideptsX;
    %insidemapY(find(insideptsY)) = i;
    %insidemapX(insidepts) = i;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [Yout,VYout] = alignImages(Y,VX,VY)

% Aligns, rotates and scales Y so VY matches VX as best as possible.
%

fprintf('Aligning images...')

mvX = mean(VX);
mvY = mean(VY);

szY = [size(Y,1), size(Y,2)];

% DV = mvX - mvY;
cVx = (VX - mvX(ones(length(VX),1),:));
cVy = (VY - mvY(ones(length(VY),1),:));

A = cVx'*cVy*(cVy'*cVy)^-1;

[U,V,W] = svd(A);

%Eliminate skew from A
V = eye(size(V))*diag(mean(diag(V)));
A = U*V*W;


[pixmapI,pixmapJ] = ind2sub(szY,1:prod(szY));

mapto = A*([pixmapJ;pixmapI] - mvY(ones(prod(szY),1),:)') + mvX(ones(prod(szY),1),:)';
VYout = cVy*A' + mvX(ones(size(VY,1),1),:);

[mshj,mshi] = meshgrid(1:szY(2),1:szY(1)); 
 
for i = 1:size(Y,3)
    
    y = Y(:,:,i);
    interpfun = TriScatteredInterp(mapto(1,:)',mapto(2,:)',y(:));
    Yout(:,:,i)=reshape(interpfun(mshj(:),mshi(:)),szY(1:2));

%     Yout(:,:,i) = griddata(reshape(mapto(1,:),szY),reshape(mapto(2,:),szY),y,mshj,mshi,'nearest',{'QJ'});

end

Yout(Yout > 1) = 1; Yout(Yout < 0) = 0;


