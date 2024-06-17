function pt = tri2patch(vert,tri,C, varargin)

%plots a triangular mesh as a patch
vx = vert(:,1);
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/tri2patch.m $
% $Revision: 37 $
% $Date: 2011-06-04 23:17:29 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

vy = vert(:,2);
vz = vert(:,3);
if nargin < 3 
    vz =vert(:,3);
    C = vz(tri);
elseif length(C) == size(vert,1);    
    C = mean(C(tri),2);
end


pt = patch(vx(tri)',vy(tri)',vz(tri)',C',varargin{:});


