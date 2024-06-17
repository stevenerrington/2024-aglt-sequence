
function  morphim = alignImages(morphim,align)

% Aligns, rotates and scales Y so VY matches VX as best as possible.
%


if nargin >1
    VX=align;
else
    VX = morphim.images(1).vert;
end

nim = length(morphim.images);
for i = 2:nim

    Y = morphim.images(i).X;
    VY = morphim.images(i).vert;

    fprintf('\nAligning image %i to 1...',i)

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

    morphim.images(i).vert = VYout;
    morphim.images(i).X = Yout;
end