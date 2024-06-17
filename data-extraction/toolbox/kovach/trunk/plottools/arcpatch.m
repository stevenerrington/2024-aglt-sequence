
function p = arcpatch(pos,theta,r)

% function p = arcpatch(pos,theta,r)
% Draws arcs centered at pos subtending theta(1) to theta(2) and r(1) to
% r(2).

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/plottools/arcpatch.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

if min(size(pos)) == 1
    pos = pos(:);
    theta = theta(:);
    r = r(:);
end

angsamp = 1;

cols = ['b','r','g','c','m','y'];
if theta(2) < theta(1)
    theta(2) = 2*pi-theta(2);
end

for i = 1:size(pos,2)
    
    th = (theta(1,i):angsamp:theta(2,i))'./180*pi;
    
    cs = [cos(th),sin(th)];
    
    vertx = cat(1,cs*r(1,i), flipud(cs)*r(2,i)) + repmat(pos(:,i)',2*length(th),1);
    
    p(i) = patch(vertx(:,1),vertx(:,2),cols(mod(i-1,length(cols))+1));
    
end

