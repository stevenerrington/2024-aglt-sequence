% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/makemesh.m $
% $Revision: 396 $
% $Date: 2013-10-28 10:26:37 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------


function [vert,tri] = makemesh(T,n,rect)

if nargin < 3
    rect = false;
end


rootv = [0 0 ; 1 0 ; 0 1 ; 1 1 ];
vert = kron(ones(n^2,1),rootv) + kron(ones(n,1),kron((1:n)'-1,repmat([1 0 ],4,1))) + kron(kron((1:n)'-1,ones(4,1)),repmat([0 1 ],n,1));
rootri = [1 2 3; 2 4 3];
tri = repmat(rootri,n^2,1) + 4*kron((1:n^2)'-1,ones(2,3));


if ~rect
    discv = vert(:,1) > n-vert(:,2);
    disctri = sum(discv(tri),2)>0;
    tri = tri(~disctri,:);
end


[unqv,unqind,memind] = unique(vert(tri(:),:),'rows');
vert = unqv./n;
tri = reshape(memind,size(memind,1)/3,3);





uu = [0 0  1; 1 0  1; 0 1  1]';

A = T/uu;

vert(:,end+1)= 1;

vert = vert*A';
