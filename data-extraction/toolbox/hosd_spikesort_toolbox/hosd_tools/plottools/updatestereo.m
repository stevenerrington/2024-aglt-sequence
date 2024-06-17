

function updatestereo(ax1,ax2,sep)


%
% update3d(ax1,ax2,sep)
%
% Positions the cameras in two axes, ax1 and ax2 to give a stereoscopic
% view (assuming the plots are the same), with eyes seperated by a given
% proportion of the distance to the midpoint (default = 0.05).
%
% See also PLOTSTEREO

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


if nargin < 3
    sep = .05;
end

axs = [ax1 ax2];


% set(axs,'DataAspectRatioMode','manual');
% set(axs,'PlotBoxAspectRatioMode','manual');
% set(axs,'CameraViewAngleMode','manual');


set(ax2,'cameraviewangle',get(ax1,'cameraviewangle'))
set(ax2,'dataaspectratio',get(ax1,'dataaspectratio'))
set([ax1,ax2],'dataaspectratiomode','manual','plotboxaspectratiomode','manual','cameraviewanglemode','manual')

cra = ismember(axs,gca);
cra = cra | all(~cra);
cp = get(axs,'CameraPosition');
cp = cra*cat(1,cp{:})./sum(cra);
ct = get(axs,'CameraTarget');
ct = cra*cat(1,ct{:})./sum(cra);
cu = get(axs,'CameraUpVector');
cu = cra*cat(1,cu{:})./sum(cra);
cva = get(axs,'CameraViewAngle');
cva = cra*cat(1,cva{:})./sum(cra);


set(ax1,'PlotBoxAspectRatio',[1 1 1])
set(ax2,'PlotBoxAspectRatio',[1 1 1])
axis manual
% axw = diff(axis);
% axw = axw(1:2:end);
% axw axdi,
rang = cross(cu,cp-ct);

cpnew = ones(length(axs),1)*cp+ ([-1 ones(1,length(ax2))].*~cra)'/sum(cra)*rang*sep;
% ctnew = [1 1]'*ct+ ([-1 1].*~cra)'/sum(cra)*rang*sep;

set(ax1,'cameraposition',cpnew(1,:),'cameraviewangle',cva);
set(ax2,'cameraposition',cpnew(2,:),'cameraviewangle',cva);
%axis(axs(cra),axis(axs(~cra)))
axis(axs(~cra),axis(axs(cra)))
% set(ax1,'cameraTarget',ctnew(1,:));
% set(ax2,'cameraTarget',ctnew(2,:));

