
function [M,MedVs] = makepic( morphim, pcnts , varargin ) 

% Function to generate morphed image. morphim is the morph image structure
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



XX = cat(4,morphim.images.X);
VX = cat(3,morphim.images.vert);
% YY = Mstruc.Y;
% VX = Mstruc.VX;
% VY = Mstruc.VY;

nim = length(morphim.images);

if ~isscalar(align)
    morphim = alignImages(align);    
elseif align
    morphim= alignImages;
end

if length(pcnts) < nim;
    pcnts = [1-sum(pcnts), pcnts];
end
szeX = size(XX);

[pixmapI,pixmapJ] = ind2sub(szeX(1:2),1:prod(szeX(1:2)));

mapto = cat(3,pixmapI,pixmapJ);

for i = 2:nim

    if size(morphim.images(i).mapto,1) == 1
      [I,J] = ind2sub(imsize(morphim.images(i)),morphim.images(i).mapto);
      mapto(i,:,:) = [I;J]';
    else
      mapto(i,:,:) = morphim.images(i).mapto';
    end

    

%     if any(any([morphim.images.detA] < 0) )&& block_neg_det
%     %    [xx,tind] = morphim.tessmap(Mstruc,pcnts + randn*.00001);  %%% Joggle slightly 
%        [xx,tind] = morphim.tessmap(pcnts); 
% 
%        nd = any(sum(tind,2)>1,2);
%        if any(nd)
%            if length(shadingscale)==1
%                shadingscale=shadingscale*ones(szeX);
%            end
%            shadingscale(nd) = foreground-1;
% 
%              shadingscale = ones(2,1)* shadingscale(:)';
% 
%        end
% 
%     end
end

if gray 
    XX = mean(XX,3);
end
    
% pixmap = [pixmapI;pixmapJ];


szmsk = size(mask);

if szmsk(3) > 1 || szmsk(2) > nim
    mask = reshape(mask,prod(szmsk(1:2)),szmsk(3));
end

if size(mask,2) < nim
%     mask = ones(2,1)*mask(:)';
    mask = [1-sum(mask,2),mask];
end    

M = {};
MedVs = {};
j = 0;
fprintf('\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t')
for j = 1:size(pcnts,1)
   
%    j = j+1;
   
   s = pcnts(j,:);
   
   fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') 
   fprintf('%2.0d of %2.0d ( %s morphing )',j,length(pcnts),sprintf('%3.0f%% ',s*100))
  
   MedVs = VX(:,:,1)*(1-sum(pcnts));
   for i = 2:nim
       
       MedVs = MedVs + pcnts(i-1)*VX(:,:,i);
   end
%    MedVs{end+1} = VX + s*(VY - VX);
   for cl = 1:size(XX,3)


%         sizem = ceil(szeX(1:2) + s * (szeY(1:2) - szeX(1:2)));

        X = reshape((XX(:,:,cl,:)),prod(szeX(1:2)),nim);

    %     [pixmapMI,pixmapMJ] = ind2sub(sizem,1:prod(sizem));

       
        if ~isequal(warpingscale , .5)
            warps = warpingscale;
            sz=size(warps);
            if size(warps,3) >1 || size(warps,2) > nim
                warps = reshape(warps,prod(sz(1:2)),sz(3));
            end
            if size(warps,2) < nim
                warps = [1-sum(warps,2),warps];
            end    

%         elseif warpingscale ~=.5
%             warps = 1;
        else
             warps = s;
        end

        if ~isequal(shadingscale ,.5)
            shades = shadingscale;
              sz=size(shades);
            if size(shades,3) >1 || size(shades,2) > nim
                shades = reshape(shades,prod(sz(1:2)),sz(3));
            end
            if size(shades,2) < nim
                shades = [1-sum(shades,2),shades];
            end   
%         elseif shadingscale == 1
%             shades = 1;
        else
             shades = s;
        end

        
           wghtsP = mask*diag(warps);
           wghtsP = wghtsP./sum(wghtsP,2);
           
           wghtsZ = mask*diag(shades);
           wghtsZ = wghtsZ./sum(wgthsZ);
            
          MedPs(1,:) =  wghtsP*mapto(:,:,1);
          MedPs(2,:) =  wghtsP*mapto(:,:,2);
          
          maptolin = sub2ind(szeX(1:2),mapto(:,:,1),mapto(:,:,2));
          MedZs(1,:) =  wghtsZ*mapto(:,:,1);
          MedZs(2,:) =  wghtsZ*mapto(:,:,2);
          
          
%           MedPs = ( (1 - mask )*( 1 - warps ) * ( 1 - warpingscale ) .* [pixmapI;pixmapJ] + mask * warps * warpingscale .* mapto)...
%                   ./ ( mask .* warps * warpingscale + (1 - mask).*(1 - warps) * (1 - warpingscale) );
% 
%           MedZ =  ( (1 - mask(1,:)' )*( 1 - shades ) * ( 1 - shadingscale(1,:)' ) .* X(:) + mask(1,:)' * shadingscale(1,:)' * shades .* Y(sub2ind(szeY(1:2),round(mapto(1,:)),round(mapto(2,:)))'))...
%                   ./ (mask(1,:)' .* shades * shadingscale(1,:)' + (1 - mask(1,:)').*(1 - shades) * (1 - shadingscale(1,:)') );

        if OthogonalWarping    %Warping orthogonal to difference vector of pixel positions
%             orthvec = [0 -1 ; 1 0]*(mapto-[pixmapI;pixmapJ]) *orthscale;
%             MedPs = MedPs + orthvec;                
            warning('Orthogonal waping is not currently implemented.')
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

          badpix = isnan(mz);
          mz(badpix)=0;
%           if ~isempty(badpix)
%               mz(badpix) = 0;
%           end

%         if ~isempty(badpix)
% 
%     %         [bdi,bdj] = ind2sub(sizem,badpix);
%     %         
%     %         for bdind = 1:length(bdi) 
%     %             bdis = max(bdi(bdind)-1),1):min(bdi(bdind)+1,sizem(1)); 
%     %             bdjs = max(bdj(bdind)-1),1):min(bdj(bdind)+1,sizem(1)); 
%     %             bdpix = meshgrid(
%     %          try
%     %             pixmeans = diag(mz(bdi+1,bdj) + mz(bdi-1,bdj) + mz(bdi,bdj-1) + mz(bdi,bdj+1) + mz(bdi-1,bdj-1) + mz(bdi+1,bdj+1)...
%     %                  + mz(bdi+1,bdj-1) + mz(bdi-1,bdj+1))/8;
%     %             mz(sub2ind(sizem,badpix)) = pixmeans;
%     %         end        
% 
%             %badpix = find(isnan(mz));
%             mz(badpix) = (( 1 - shades ) * ( 1 - shadingscale(1) ) * X(badpix) + shadingscale(1) * shades * Y(badpix))./(shades * shadingscale (1)+ (1 - shades) * (1 - shadingscale(1)) );
% 
%         end  


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






