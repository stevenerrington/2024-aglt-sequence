%%% One-liner utilities
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/one_liners.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

% matrix plot functions: avoid the annoyance of giving each column as a
% separate argument.

mtrisurf = @(tri,X,varargin)  trisurf(tri,X(:,1),X(:,2),X(:,3),varargin{:});
mtrimesh = @(tri,X,varargin)  trimesh(tri,X(:,1),X(:,2),X(:,3),varargin{:});
trsurf = @(tr, varargin)mtrisurf(tr.Triangulation,tr.X,varargin{:});
trmesh = @(tr, varargin)mtrimesh(tr.Triangulation,tr.X,varargin{:});
%Cross product
crpr = @(x,tr)cross(x(tr(:,2),:)-x(tr(:,1),:),x(tr(:,3),:)-x(tr(:,1),:));
%compute area of a facet.
facearea = @(tr)sqrt(sum(crpr(tr.X,tr.Triangulation).^2,2));
mplot3 = @(x,varargin)plot3(x(:,1),x(:,2),x(:,3),varargin{:});
m2plot3 = @(x,y,varargin)  plot3(cat(2,x(:,1),y(:,1))',cat(2,x(:,2),y(:,2))',cat(2,x(:,3),y(:,3))',varargin{:});
mplot = @(x,varargin)plot(x(:,1),x(:,2),varargin{:});
mquiver = @(X,V,varargin) quiver3(X(:,1),X(:,2),X(:,3),V(:,1),V(:,2),V(:,3),varargin{:});
shrink = @(x,a) a*(x-repmat(mean(x),size(x,1),1))+repmat(mean(x),size(x,1),1);
trshrink = @(tr,a)TriRep(tr.Triangulation, shrink(tr.X,a));
wire = @(h)set(h,'edgealpha',1,'facealpha',0);
face = @(h)set(h,'edgealpha',0,'facealpha',1);
pad = @(x) cat(2,x,ones(size(x,1),1));