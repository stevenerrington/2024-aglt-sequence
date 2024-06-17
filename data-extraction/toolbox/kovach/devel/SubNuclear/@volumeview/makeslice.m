function [SL,edge,trmat,hn] = makeslice(me,nv,pt,rot)

%%% return a slice normal to nv containing pt rotated about nv by rot


% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/@volumeview/makeslice.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 4
    rot = 0;
end

nv = nv./norm(nv);

fig = figure;
set(fig,'visible','off')

sz = size(me.Vol);
sz = sz([2 1 3]);
h = pcolor((0:sz(2)-1)+2,(0:sz(1)-1)+2,zeros(size(me.Vol(:,:,1)))');
set(h,'zdata',ones(size(me.Vol(:,:,1)))*round(pt(3)))

rax = cross( nv,[0 0 1]);
if norm(rax) == 0
    rax = [0 0 1];
    th = 0;
    M = eye(3);
else
    rax = rax./norm(rax);
    th = acos([0 0 1]*nv');
    plax = cross([0 0 1],rax); %axis in plane before rotation
    plax2 = cross(nv,rax);  % axis in plane before rotation
    M = [0 0 1; plax; rax]\[nv;plax2;rax]; % Rotation matrix
end
% axv = [1 0 0; 0 1 0]*M; % axis vectors for the new plane
 rotate(h,rax,-th/pi*180,pt);
if rot > 0
    rotate(h,nv,-rot/pi*180,pt)
end

 xd = get(h,'xdata');
 yd = get(h,'ydata');
 zd = get(h,'zdata');


hn = slice(me.Vol,yd,xd,zd);

% xy = [xd([1 end]); yd([1 end]); zd([1 end])]';


SL = get(hn,'cdata')';

x0 = [0 0  ;  0 sz(2)-1 ;sz(1)-1 0 ; sz(1:2)-1 ; 0 0 ]+2;
ln = @(x)x(:);
corn = @(x)ln(x([1 end],[1 end]));

x1 = [corn(xd),corn(yd), corn(zd), ones(4,1)];
x1 = cat(1,x1,x1(1,:) + [nv,0]);
trmat = x1\x0;

if nargout <4
    delete(fig)
else
    set(fig,'visible','on')
end
edge = x1;
