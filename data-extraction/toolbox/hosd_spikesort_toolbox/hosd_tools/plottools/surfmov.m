function avi = surfmov(X,Y,Z,varargin)

%   avi = surfmov(X,Y,Z,KEYWORDS)
%   SURFMOV creates a movie with a pcolor plot of X aligned with Y, and
%   superimposed on Z. Each frame corresponds to a row of X.
%
%   Keywords:
%
%   'Filename':       Name of file to be stored. If an identically named
%                     file exists, the movie will be stored as FILENAME_N.AVI,
%                     where N is the smallest integer for which a file with
%                     the same name doesn't already exist. If the filename
%                     isn't specified the movie will be stored as
%                     brainmovie_N.avi. 
%   'Transparency':   transparency of the plot. 0 to 1.
%   'GridSize':       followed by [M N]. Maps each row of X onto an M by N grid.
%                     Default is 8 by 8. N is the number of contacts in a row and M 
%                     the number of rows, where "row" means a line of contiguously 
%                     numbered contacts. 
%   'KeepAll':        Retain axis positioning information from previous
%                     call to SUPERAX.
%   'Size':           Factor by which figure size is reduced from full
%                       screen when capturing a frame. Smaller size ->
%                       smaller file size.
%   'Colorbar':       Adds colorbar at edge of frame.
%   'Caxis':          Use given color axis.
%   'StimOnset':      Onset of Stimulus in Msec.
%   'SamplingFreq':   sampling frequency.
%   'Start':          Time in msec at which to begin recording relative to onset of stimulus in msec.
%   'Finish':         Time at which to end recording
%   'Jump':           Use every ith data time point per movie frame. Default is 3.
%   'Shading':        Default is 'interp'. 'Flat' or 'faceted' are also valid options (see HELP SHADING).
%   'SamplingFreq':   Sampling frequency
%   'Compression':    Video compression

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C. Kovach 2002
%Modified 9/15/02

persistent old_ax2

%-------movie settings----
compression = 'indeo5';

%------Default Keyword values---------
alpha = .6;
varg = varargin;
grid_size = [8 8];
pass = {};
size_reduction = .6;
cbar = 0;
fname = [];
cax = [.8*min(min(X)), .8*max(max(X))];
stonset = 0;
Fs = 2000;
start = [];
fin = [];
skip = 3;
notext = 0;
backgr = 1;
%---------------------------------
i = 1;
while i <= length(varg)
  switch lower(varg{i})
      case 'transparency'
          alpha = 1 - varg{i+1};
          i = i+1;
      case 'gridsize'
          grid_size = varg{i+1};
          i = i+1;
      case 'keepall'
          pass{end + 1} = 'keepall';
      case 'size'
          size_reduction = varg{i+1};
           i = i+1;
       case 'shading'
           pass(end + 1:end+2) = {'shading', varg{i+1}};
           i = i+1;
      case 'colorbar'
          cbar = 1;
      case 'filename'
          fname = varg{i+1};
          i = i+1;
      case 'caxis'
          cax = varg{i+1};
          i = i+1;
      case 'stimonset'
          stonset = varg{i+1};
          i = i+1;
      case 'samplingfreq'
          Fs = varg{i+1};
          i = i+1;
      case 'start'
          start = varg{i+1};
          i = i+1;
      case 'finish' 
          fin = varg{i+1};          
        i = i+1; 
      case 'jump'
          skip = varg{i+1};
          i = i+1;
      case 'supresstext'
          notext = 1;
      case 'nobackgr'
          backgr = 0;
      case 'compression'
          compression= varg{i+1};          
        i = i+1;
      otherwise
          error(['''',varg{i},''' isn''t a valid keyword'])
  end
  i = i+1;
end

pass(end + 1:end+2) = {'size', size_reduction};

if isempty(fname)
    fname = 'brainmovie';  
end


%Creates a new filename if an old one already exists
j = [];
addon = [];
k = 1;

while k
    nuname = [fname,addon];
    k = (exist([nuname,'.avi'],'file') ~= 0);
    if isempty(j)
        j = 0;
    end
    j = j + 1;
    addon = ['_',num2str(j)];
end

%reuses previously used axes
if ~isempty(old_ax2)
    ax2 = old_ax2;
else
    ax2 = axes;
end 

 %shows background
if backgr
 s = imagesc(Z); colormap gray, axis image, axis off
 if size(Z,3) < 3
     Z = cscale(s,'rgb');
 end


 %Shows image with which foreground will be aligned

 s = imagesc(Y); colormap gray, axis image, axis off
 if size(Y,3) < 3
     Y = cscale(s,'rgb');
     imagesc(Y);
 end
 

 %Generates a test-grid 
 test_grid = zeros(size(X(1,:)));
 test_grid(1) = 1;
 test_grid(end) = -1;
 test_grid(grid_size(2)) = .5;

 c = square(test_grid,grid_size(1), grid_size(2));


 %obtains positioning data with superax.
 axis off
 [ax1,ax2,reax,review, refl] = superax(c,ax2,pass{:});
 s = get(ax1,'children');

 %pause       %pause to check if good superposition
end

axes(ax1)
caxis(cax)
if cbar
 cb = colorbar;
 cb_pos = get(cb,'position');
end

if backgr
set(ax1,'position',reax,'view',review)
end


%delete(ax1)

%Opens avi-file
avi = avifile(nuname,'compression',compression);

%Shows background on axis 2
if backgr
 pos = get(ax2,'position');
 axes(ax2), cla
 imagesc(Z), axis image, axis off
end

%makes and positions colorbar

caxis(cax)
if cbar
 cb = colorbar;
 cb_child = get(cb,'children');
 set(ax2,'position',pos)
 set(cb,'position',[pos(1)+pos(3)-5*cb_pos(3),pos(2) + .05*cb_pos(4),cb_pos(3:4)])
 set(cb_child,'alphadata',alpha)
end

%Generates and positions title

 title_handle = title(''); 
 set(title_handle,'interpreter','none','units','normalized','fontsize',15)
 pos = get(ax2,'position');
 title_position = get(title_handle,'position');
 set(title_handle,'position',[title_position(1)- .2*pos(3), pos(2) - .1*pos(4) ,title_position(3)])


%Changes start-time in msec to frame number
if isempty(start)
   framestart = 1;
else
    framestart = round((start+stonset)*Fs/1000);
end

%Same for end-time
if isempty(fin)
   framefin = size(X,1);
else
    framefin = round((fin+stonset)*Fs/1000);
end

%makes movie
for t = framestart:skip:framefin;
    
    
    c = square(X(t,:),grid_size(1),grid_size(2),s, refl);
    
    set(s,'CData',c);
    caxis(cax);
   if ~notext   
    set(title_handle,'string',['           ',fname,...
        '                           ',...  
        'Time: ',sprintf('%5.1f',t*1000/Fs - stonset),'ms post-stimulus'],'fontsize',10)
    
   elseif ~ischar(stonset)
        set(title_handle,'string','Time: ',sprintf('%5.1f',t*1000/Fs - stonset,'ms post-stimulus'),'fontsize',10)
   end
    axes(ax1)
    if cbar
     axes(cb)
    end
    frame = getframe(ax2);
    avi=addframe(avi,frame);
    
end

delete(gcf)

if nargout == 0
    avi = close(avi);
end
%-----------------------------------
function c = square(x,m,n,s, refl);

if nargin < 4
    refl = 0;
end


%Converts x into a rectangle of dimensions nXm.
c = reshape(x',n,m)';

if refl
    c = fliplr(c);
end

if nargin > 3
   c= flipud(c);           %Flips ud inorder to make pcolor image with same orientation as imagesc
   cd = get(s,'cdata');

   cd(1:size(c,1),1:size(c,2)) = c;  %intended to preserve padding if shading is other than interp

   c = cd;
end




