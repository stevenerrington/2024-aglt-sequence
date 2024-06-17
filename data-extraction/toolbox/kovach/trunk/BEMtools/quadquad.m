function I = quadquad(bem,varargin)

% Quadrature (numerical integration) of a quadratic approximation on a tesselated surface

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/quadquad.m $
% $Revision: 37 $
% $Date: 2011-06-04 23:17:29 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

if isnumeric(varargin{1})
    Xvert = varargin{1};
    Xedge = varargin{2};
    
else
    intfun = varargin{1};
    Xvert = intfun(bem.pnt);
    Xedge = intfun((bem.pnt(bem.edges(:,1))+bem.pnt(bem.edges(:,2)))/2);
end


N = tessquadinterp(bem,Xvert,Xedge);

pnt = bem.pnt;
%Edge vectors
E1 = pnt(bem.tri(:,2),:) - pnt(bem.tri(:,1),:);
E2 = pnt(bem.tri(:,3),:) - pnt(bem.tri(:,1),:);

normvec = cross(E1,E2);
area = sqrt(sum(normvec.^2,2))./scale;

I = area.*sum(N,2)./6;







