function [IN,ON,PP] = intersects(X1,X2,tri,pnt)

%Finds tesselations intersected by lines connecting the rows in X1 and X2.
%
%  I(i,j) is 1 if the line connecting X1(j,:) to X2(j,:) intesects facet i

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/intersects.m $
% $Revision: 191 $
% $Date: 2013-04-24 09:53:04 -0500 (Wed, 24 Apr 2013) $
% $Author: ckovach $
% ------------------------------------------------


%normal vectors
normvec = (cross(pnt(tri(:,2),:) - pnt(tri(:,1),:),pnt(tri(:,3),:)- pnt(tri(:,1),:)))';
area = sqrt(sum(normvec.^2));
unorm = normvec./repmat(area,size(normvec,1),1);
unorm(isnan(unorm))=0;
% cent = (pnt(tri(1,:),:) + pnt(tri(2,:),:)+pnt(tri(3,:),:))./3;
base = pnt(tri(:,1),:)';
PL1 = pnt(tri(:,2),:)' - base;
PL2 = pnt(tri(:,3),:)' - base;
% INV1 = PL2*[0 -1; 1 0]./area;
% INV2 = PL2*[0 1; -1 0]./area;
SPL = sparseblock(cat(2,PL1(:),PL2(:),unorm(:)),3,'transpose');
SPL = SPL + speye(length(SPL))*eps;
nvert = size(tri,1);

IN = false(size(tri,1),size(X1,1));
ON = false(size(tri,1),size(X1,1));
if nargout > 2
    PP = zeros([size(tri,1),size(X1,1) 3]);    
end

for i =1:size(X1,1)
    
    DX = X2(i,:)-X1(i,:);
    %project Xs onto normal dimension for each tesselation
    xn = X1(i,:);
    
    p1 = xn*unorm - sum(unorm.*base);
    dp = DX*unorm +eps;

    s = p1./dp; % 

    PR = (repmat(xn',1,nvert) - DX'*s) -  base;

    %project into unit spaces
    U = reshape((SPL + 1e-6*speye(size(SPL)))\PR(:),3,nvert);

% Test
    IN(:,i) = (sign(p1) ~= sign(p1+dp)) & U(1,:)>=0 & U(2,:)>=0& (U(1,:) <= 1-U(2,:));
    ON(:,i) = (p1==0) & U(1,:)>=0 & U(2,:)>=0& (U(1,:) <= 1-U(2,:)) + 2*(p1+dp==0) & U(1,:)>=0 & U(2,:)>=0& (U(1,:) <= 1-U(2,:));
    
    if nargout>2
        PP(IN(:,i),i,:) = PR(:,IN(:,i))' + base(:,IN(:,i))';
    end
end

