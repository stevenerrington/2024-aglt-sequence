
function tmp = template206


%Create templates for 180

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/template206.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------


z = 0;

dx = 5;
[Y,X] = meshgrid((1:12)+10,(1:8)-10); %research grid
XY = [X(:),Y(:)]*dx;
Ds{1} = dfun(XY,XY) ;
XYs{1}(:,[1 3]) = XY;
XYs{end}(:,2) = z;


dx = 10;
[Y,X] = meshgrid((1:8) +10 ,(1:4)-10); %clinical grid
XY = [X(:),Y(:)]*dx;
Ds{2} = dfun(XY,XY);
XYs{2}(:,[1 3]) = XY;
XYs{end}(:,2) = z;


Ds{1} = blkdiag(Ds{1:2});
Ds{1}(Ds{1}==0) = -1;
Ds{1} = Ds{1}.*(1-eye(length(Ds{1})));
XYs{1} = cat(1,XYs{1:2});
Ds(2) = [];
XYs(2) = [];

nstrip = 3;
for k = 1:nstrip
    dx = 10;
    [Y,X] = meshgrid(2*k,(1:4) ); %4-contact strip
    XY = [X(:),Y(:)]*dx;
    Ds{end+1} = dfun(XY,XY);
    XYs{end+1}(:,[1 3]) = XY;
    XYs{end}(:,2) = z;

end

n16grid = 2;
for k = 1:n16grid
    dx = 5;
    [Y,X] = meshgrid(2*nstrip+4*k+[1:2],(1:8)); %16-contact grid
    XY = [X(:),Y(:)]*dx;
    Ds{end+1} = dfun(XY,XY);
    XYs{end+1}(:,[1 3]) = XY;
    XYs{end}(:,2) = z;

end



n10strip = 1;
for k = 1:n10strip
    dx = 10;
    [Y,X] = meshgrid(2*nstrip+4*n16grid,(1:10)); 
    XY = [X(:),Y(:)]*dx;
    Ds{end+1} = dfun(XY,XY);
    XYs{end+1}(:,[1 3]) = XY;
 XYs{end}(:,2) = z;
end


ntpgrid1 = 1;
for k = 1:ntpgrid1
    dx = 10;
    [Y,X] = meshgrid(1:6,(1:6)); %temporal pole will be modeled as separate strips
    keep = (X<2|Y<6) & (X<3|Y<5); 
    
    Y = Y+2*nstrip+4*n16grid + 6*n10strip;
    XY = [X(keep),Y(keep)]*dx;
    Ds{end+1} = dfun(XY,XY);
    XYs{end+1}(:,[1 3]) = XY;
    XYs{end}(:,2) = z;
end
%Template distances
T = blkdiag(Ds{:});
T(T==0) = -1;
T = T.*(1-eye(length(T)));
T(end+1,:) = -1;
T(:,end+1) = -1;

tmp.Ts = Ds;
tmp.T = T;
tmp.TXYs = XYs;

