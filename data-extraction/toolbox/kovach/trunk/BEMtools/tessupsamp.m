function [Y,X,tri] = tessupsamp(bem,N,n)

%Returns quadratic interpolation upsampled by a factor of n for each tesselation 

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/tessupsamp.m $
% $Revision: 37 $
% $Date: 2011-06-04 23:17:29 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 3
    n = 4;
end

[w,wtri] = makemesh([0 0 0;1 0 0;  0 1 0]',n);

%linear basis functions
PHI =  [1-w(:,1)-w(:,2),w(:,1),w(:,2)];

%quadratic basis
PHQ = kron(PHI,ones(1,3)).*kron(ones(1,3),PHI);
PHQ = PHQ(:,[1 5 9 2 6 3]);

Y = PHQ*N';

Y = Y(:);

%edges
base = bem.pnt(bem.tri(:,1),:);
E1 = bem.pnt(bem.tri(:,2),:)-base;
E2 = bem.pnt(bem.tri(:,3),:)-base;

X = kron(E1,w(:,1)) + kron(E2,w(:,2)) + kron(base,ones(size(w,1),1));

tri = repmat(wtri,size(E1,1),1) + kron((1:size(E1,1))'-1,ones(size(wtri,1),3))*size(w,1);




