
function varargout = plotstereo(x,y,z,varargin)


% plotstereo(x,y,z,varargin)
%
% Like plot3, except it generates a stereoscopic view in two subplots.
%
%
% See also SETSTEREO and UPDATESTERO

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


sep = .05;
anaglyph = false;
anaglyphcol = [1 0 0; 0 0 1];
i = 1;
while i <= length(varargin)
    
    switch lower(varargin{i})             
        case 'sep' % specify separations
            sep  = varargin{i+1};
            i = i+1; 
        case 'anaglyph' % specify separations
            anaglyph  = varargin{i+1};
            i = i+1; 
    otherwise
           
            break
%             error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;    
end

plotargs = varargin(i:end);

if anaglyph
    sbpl1 = {1,1,1};
    sbpl2 = {1,1,1};
    plarg1 = [plotargs,{'color',anaglyphcol(1,:)}];
    plarg2 = [plotargs,{'color',anaglyphcol(2,:)}];
    hold on
else
    sbpl1 = {1,2,1};
    sbpl2 = {1,2,2};
    plarg1 = plotargs;
    plarg2 = plotargs;
end

ax1 = subplot(sbpl1{:});

pl = plot3(x,y,z,plarg1{:});

ax2 = subplot(sbpl2{:});

pl = [pl;plot3(x,y,z,plarg2{:})];


set(ax1,'PlotBoxAspectRatio',[1 1 1])
set(ax2,'PlotBoxAspectRatio',[1 1 1])


setstereo(ax1,ax2,sep)

axis(ax1,'vis3d')
axis(ax2,'vis3d')


if nargout > 0
    varargout{1} = pl;
end