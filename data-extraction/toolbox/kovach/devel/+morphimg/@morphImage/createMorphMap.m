function createMorphMap(morphim,varargin)

%
% Constructor for morphImage class
%

import morphimg.*

morphim.images = imagedat;
morphim.images(1) = []; 

if isequal(class(varargin{1}),'morphImage') 
    
    morphim= varargin{1};

else
    for k = 1:length(varargin)
        if ischar(varargin{k})
            varargin = varargin(k:end);
            break            
        elseif isequal(class(varargin{1}),'imagedat') 
            morphim.images(k)= varargin{k};
        else
            morphim.images(k) = imagedat(varargin{k});
        end
    end
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
            for k = 1:length(morphim.images)
                if ischar(varargin{k+i})
                    break
                    i=i-1;
                else
                    morphim.images(k).vert = varargin{k+i};                    
                end                
            end
            i=i+k+1;
            
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

% tempsize = [max([size(X,1),size(Y,1)]),max([size(X,2),size(Y,2)]),max([size(X,3),size(Y,3)])];
tempsize = [max(morphim.images.imsize),3*max(morphim.images.truecolor)];
template = zeros(tempsize);

for i = 1:length(morphim.images)
    if ~isequal(morphim.images(i).imsize,tempsize)
       T = template;
       X=morphim.images(i).X;
       T(1:size(X,1),1:size(X,2),:) = X;
       morphim.images(i).X = T;
    end
    clear T
end



nv = morphim.images.nvert;
mxnv = max(nv);

if length(Tesswgt) < nv
    Tesswgt =[Tesswgt (1-Tesswgt)];
end

% fmx = find(nv==mxnv);fmx(2:end) =[];

strs = {' has','s have'};

for i =1:length(nv)
    if nv(i)<mxnv
        str = strs{(mxnv-nv(i) > 1)+1};
    
        warning('Number of Vertex points doesn''t match. %i point%s been added to image %i',mxnv-nv(i),str)
        morphim.images(i).vert(end+1:mxnv,:) = morphim.images(mxnv).vert(size(VYin,1)+1:end,:);
    end
    
end    


if getpoints
    morphim.getvertices(fixax,Tesswgt);
end

if ~isscalar(align)
    morphim.alignImages();    
elseif align
     morphim.alignImages();
end

fprintf('Computing pixel map... |')

VTess = morphim.images.weightedsum( Tesswgt ) ;

Tess = delaunay(VTess(:,1),VTess(:,2)); % delaunay tesselation of fiducial points;

szeX = morphim.images(1).imsize;
[pixmapI,pixmapJ] = ind2sub(szeX(1:2),1:prod(szeX(1:2)));
pixmap = [pixmapI;pixmapJ];
% [pixmapYI,pixmapYJ] = ind2sub(szeY(1:2),1:prod(szeY(1:2)));

mapto = zeros(2,size(pixmapI,2));
%insidemapY = zeros(szeY(1:2));
insidemapX = zeros(szeX(1:2));
rotor = '/-\|';
detA = 0;

nim = length(morphim.images);

% VM = morphim.images.weightedsum;
VM = morphim.images(1).vert;
JI = [pixmapJ;pixmapI];

for k = 1:nim
     morphim.images(k).mapto = mapto;
end

for i = 1:size(Tess,1)
    insidepts = find(inpolygon(pixmapJ,pixmapI,VM(Tess(i,:),1),VM(Tess(i,:),2)));
    for k = 1:nim

    if k == 1
%         morphim.images(1).mapto = [pixmalJ,pixmapI]';
        morphim.images(1).mapto = [];
        morphim.images(1).detA(i,1)=1;

    else

            fprintf('\b%s',rotor(mod(i,4)+1))

            insidemapX(insidepts) = i;


            %Linear transform: y = Ax + b; Since vy = A*(vx - vx0) + vy0, and [vy1 vy2] = A*([vx1,vx2] - [vx0,vx0]) + [vy0,vy0]
            % then A = ([vy1 vy2] - [vy0,vy0]) * ([vx1,vx2] - [vx0,vx0])^-1;
            A = ([morphim.images(k).vert(Tess(i,2),:)',morphim.images(k).vert(Tess(i,3),:)'] - [morphim.images(k).vert(Tess(i,1),:)',morphim.images(k).vert(Tess(i,1),:)'])...
                * ([VM(Tess(i,2),:)',VM(Tess(i,3),:)'] - [VM(Tess(i,1),:)',VM(Tess(i,1),:)'])^-1;

            %Vertices are in XY form, hence flipud to make mapto IJ.
            morphim.images(k).mapto(:,insidepts) =  flipud(A * ([pixmapJ(insidepts);pixmapI(insidepts)] - VM(Tess(i,1),:)'*ones(1,length(insidepts)))...
                                    + morphim.images(k).vert(Tess(i,1),:)'*ones(1,length(insidepts))) -JI(:,insidepts);
            
            morphim.images(k).detA(i,1)=det(A);


    
        end
        
%         morphim.images(k).detA=detA;
        
    end
    
%     morphim.images(k).facets = facet(morphim.images(1).vert,VM,Tess);
end

morphim.Tri = Tess;

for k = 2:nim
 
% if ~eraseUnmappedPixels
 if roundmapto
      morphim.images(k).mapto = round(morphim.images(k).mapto);
 end
    
  
  sz = morphim.images(1).imsize;
  
  mapto = morphim.images(k).mapto ;
  mapto(mapto <= 0) = pixmap(mapto <= 0);
  mapto(1,(mapto(1,:) > sz(1))) = pixmap((mapto(1,:) > sz(1)));
  mapto(2,(mapto(2,:) > sz(2))) = pixmap((mapto(2,:) > sz(2)));
  
  [Iind,Jind] = ind2sub(sz,find(isnan(mapto(1,:))));
  mapto(:,isnan(mapto(1,:))) = [Iind; Jind];
  
  %else
  %unmappedX = pixmap(mapto == 0);  
  %mappedY = mapto(mapto(
  %X(unmapped(1,:),unmapped(2,:)) = 0;
  %Y(ismember(

    % morphim.X = X;
    % morphim.Y = Y;
    % morphim.VX = VX;
    % morphim.VY = VY;
    % morphim.Tess = Tess;
    % morphim.detA=detA;
    % morphim.detA=Tesswgt;
    if roundmapto
      morphim.images(k).mapto = sub2ind(size(X),mapto(1,:),mapto(2,:));
    else
      morphim.images(k).mapto = mapto-JI;
    end
end
%
% morphim.makepic = @makepicfun;
% morphim.tessmap = @tessmapfun;


