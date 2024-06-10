function [pts,transf] = fitgrid(P,MN,PX,lin,label)

% pts = fit_grid(P,[m n],PX,lin);
%
% Interpolates the points of an m x n grid based on the coordinates in matrix P.
%
% Inputs:
%
% P: Kx3 matrix
%
% [m n]: dimensions of the grid
%
% PX: Grid points corresponding to the rows of P 
%   If PX is not given or empty it defaults to
%       PX = [1 1; m 1; 1 n; m n]
%
% lin: Uses linear interpolation if lin = true and otherwise thin-plate spline
%       warping to project grid points into the space of P.
%      If lin is not given or empty, it defaults to true when k < 5 and
%      otherwise false.
%
% Outputs:
%
% pts: m*n x 3 array of interpolated points.

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/fitgrid.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

% C. Kovach 2013


MN(3) = 1;
if nargin < 3 || isempty(PX)
    PX = MN([3 3; 1 3; 3 2; 1 2]);
end

if  nargin < 4 || isempty(lin);
  lin =  size(PX,1) < 5; 
end

if nargin  <5 || isempty(label)
   label = 'Grid %i'; 
end

if isa(P,'points')
    P = cat(1,P.coord);
end
[x, y] = meshgrid(1:MN(2),1:MN(1));

if lin
    P(:,4) = 1;
    PX(:,3) = 1;
    T = PX\P;
    transf = transforms('trmat',T,'label','To Grid');
   
%     trfun = @(X) [X ones(size(X,1),1)]*T(:,1:end-1); 
    
   
else
   transf = transforms('label','To Grid');
   tpsfun = tpswarp(PX,P,1e-9); 
%    trfun = @(X) tpsfun(X);
   transf.tr = tpsfun;
   if nargout > 1
       tpsifun = tpswarp(P,PX,1e-9); 
      transf.itr = tpsifun;
   end   
end

pts = points(transf.tr([y(:) x(:) ]));
for k = 1:length(pts)
    pts(k).label = sprintf(label,k);
end
