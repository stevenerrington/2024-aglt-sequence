function N = tessquadinterp(bem,Xvert,Xedge)

%fit a quadratic function over a tesselated surface 
%
%Xvert - value at vertices.
%Xedge - value at edge midpoints.
%
% The output N, is a matrix of parameter values for the set of basis
% functions, which are the kronecker product of linear basis functions, 
% f,over each face. f has the property that fi(vj) = delta[i,j].
% 
%                     p2
%                     /\
%                    /  \
%                p4 /    \ p5
%                  /      \
%                 /________\
%               p1    p6    p3
%               
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/tessquadinterp.m $
% $Revision: 126 $
% $Date: 2012-06-06 14:31:56 -0500 (Wed, 06 Jun 2012) $
% $Author: ckovach $
% ------------------------------------------------


% C. Kovach 2009


Xedge = Xedge(:);
Xvert = Xvert(:)';

N =  Xvert(bem.tri);

for i = 1:3
    N(:,3+i)   =  4*Xedge(bem.edgemap(:,i)) - sum(Xvert(bem.edges(bem.edgemap(:,i),:)),2);
end

% bem.quadfit = N;

