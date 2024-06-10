
function [Yout,VYout,A] = rigidAlign(Y,VX,VY,noskew)

% Computes a rigid linear transformation which minimizes squared diffirence
% between VX and VY, then applies the transformation to points in matrix Y.
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/imagetools/rigidAlign.m $
% $Revision: 191 $
% $Date: 2013-04-24 09:53:04 -0500 (Wed, 24 Apr 2013) $
% $Author: ckovach $
% ------------------------------------------------


if nargin < 4
    noskew = 1;
end

fprintf('Aligning images...')

mvX = mean(VX);
mvY = mean(VY);

szY = [size(Y,1), size(Y,2)];

DV = mvX - mvY;
cVx = (VX - mvX(ones(length(VX),1),:));
cVy = (VY - mvY(ones(length(VY),1),:));

A = cVx'*cVy*(cVy'*cVy)^-1;

[U,V,W] = svd(A);

if noskew
    %Eliminate skew from A
    V = eye(size(V))*diag(mean(diag(V)));
    A = U*V*W;
end


[pixmapI,pixmapJ] = ind2sub(szY,1:prod(szY));

mapto = A*([pixmapJ;pixmapI] - mvY(ones(prod(szY),1),:)') + mvX(ones(prod(szY),1),:)';
VYout = cVy*A' + mvX(ones(size(VY,1),1),:);

[mshj,mshi] = meshgrid(1:szY(2),1:szY(1)); 
 
for i = 1:size(Y,3)
    
    y = double(Y(:,:,i));
    
%     Yout(:,:,i) = griddata(reshape(mapto(1,:),szY),reshape(mapto(2,:),szY),y,mshj,mshi,'nearest');
    interpfun = TriScatteredInterp(mapto(1,:)',mapto(2,:)',y(:),'nearest');
    Yout(:,:,i)=reshape(interpfun(mshj(:),mshi(:)),szY(1:2));

end

Yout(Yout > 1) = 1; Yout(Yout < 0) = 0;



