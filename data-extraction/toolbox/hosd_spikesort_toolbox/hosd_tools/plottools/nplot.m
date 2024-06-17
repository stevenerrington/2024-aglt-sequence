function varargout = nplot(Dat,varargin)

%axs = nplot(Dat,varargin)
%NPLOT makes a figure with subplots for each column in Dat.
%Optional Keywords:
%   'Title': title of figure.
%   'Fig': handle of figure in which to plot.
%   'ChLabels': vector of channel labels (default is the column number
%       associated with each channel in the data matrix
%   'Axes': followed by string containing command applied separately to
%       each axis. 
%   'Axis': sets axis window.
%   'Flip': flips the grid left to right. To flip up to down, use 'flip'
%           and ('rotate',180).
%   'TicksOn': followed by channel(s) for which axis tick-marks are
%       displayed.
%   'Rotate': rotates position of plots clockwise by specified degrees.
%   'White':  makes background white instead of default color (gray).
%   'Size':   [m n] dimensions of subplot (rows columns); default is 8-by-8.
%   'Superimpose': Followed by axs vector from previous call to nplot.
%                  Superimposes plots.
%   'LineStyle':  type 'help plot' for line styles
%   'LineWidth':  Width of line. 1 is default.
%   'LineColor': Defines line color with RGB value
%   'SamplingFreq': sampling frequency. Causes xticks to be in milliseconds and 'msec' to appear as xlabel. 
%   'StOnset': Creates vertical line at specified point. Units msec if
%            sampling frequency is given, otherwise units are sampling period.
%   'HorzLine': creates horizontal line at specified point of the ordinate.
%   'Flipy':  Plots positive values down and negative up.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C. Kovach 2002
varg = varargin;

i = 1;

%----------------------  Default settings
title_size = 10;
flip = 0;
theta = 0;
ax_command = {};
tickson = 57;
chlbl = [];
fig = [];
backwhite = 0;
grid_size = [8 8];
axes_handles = [];
plot_command = {};
new_plot = 1;
axlim = [];
stonset = 0;
Fs = 1;
title_height = 1;       %Height of the labels on individual sub-plots in normalized axis units.
horzline = [];
channels = [];
flipy = 0;
linecolor = [];
axis_line_style = '-';
axis_line_color = 'k';
xdata = [];
%--------------------------

if ~isempty(varargin) && isnumeric(varargin{1})
    xdata = Dat;
    Dat = varargin{1};
    Fs = [];
    i = i+1;
end

while i <= length(varg)
   switch lower(varg{i})
       case 'title'
           stim = varg{i+1};
       case 'fig'
           fig = varg{i+1};
       case 'chlabels'
           chlbl = varg{i+1};
  %    case 'labels'
  %        nchan = varg{i+1};
       case 'axes'
           ax_command = cat(2, ax_command,varg(i+1));
       case 'flip'
           flip = 1;
           i = i-1;
          
           if tickson == 57
               tickson = 64;
           end
       case 'tickson'
           tickson = varg{i+1};
       case 'rotate'
           theta = varg{i+1};
       case 'white'
           backwhite = 1;
           i = i-1;
       case 'size'
           grid_size = varg{i+1};
       case 'superimpose'
           axes_handles = varg{i+1};
           new_plot = 0;
       case 'linestyle'
           plot_command = cat(2,varg(i+1), plot_command);
       case 'linewidth'
           plot_command = cat(2,plot_command,{'linewidth'},varg(i+1));
       case 'linecolor'
           linecolor = varg{i+1};
       case 'axis'
           axlim = varg{i+1};
       case 'stonset'
           stonset = varg{i+1};
%            stonset = [stonset,varg{i+1}];
       case 'samplingfreq'
           Fs = varg{i+1};
       case 'horzline'
           horzline = [horzline,varg{i+1}];
       case 'channels'
           channels = [channels,varg{i+1}];
       case 'flipy'
           flipy = 1;
           i = i-1;
       otherwise
           error([varg{i},' is not a valid option.']);
   end         
   i = i+2;
end



%if class(Dat) == 'struct'
%    nchan = Dat.nchan;
    %elseif nargin < 5
   % nchan = size(Dat,2);
   %else 
    
   if ~isempty(channels)
       Dat = Dat(:,channels);
     if ~isempty(chlbl)
         chlbl = chlbl(:,channels);
     else
         chlbl = channels;
     end
   end
   
   if isempty(axlim) && ~isempty(Fs)
       axlim = [-stonset, size(Dat,1)./Fs-stonset,min(min(Dat)),max(max(Dat))]; 
   else

       axlim = [min(xdata),max(xdata),min(min(Dat)),max(max(Dat))]; 
   end       
   
   nchan = size(Dat,2);
   
   plotme = true(1,nchan);
   %end

%if nargin < 4
%    chlbl = 1:nchan;
%elseif isempty(chlbl)
%    chlbl = 1:nchan;
%end

%if nargin < 3
%    figure
%else
%    figure(fig)
%end

if isempty(fig) & new_plot
    fig = figure;
end

 %if class(Dat) == 'struct'
 %   eval(['G = Dat.',stim,';']);
 %else
 
 G = Dat;
    
 % end

 if size(G,2) > 200
     fprintf('/n>Sure you want to plot %i channels? (^C to exit, any key to continue)', size(G,2))
     pause
 end
 
 channum = 1:size(G,2);
 
 


 if flip == 1
     index = reshape([1:prod(grid_size)]', grid_size(2), grid_size(1))';
     index = fliplr(index);
    index = reshape(index',prod(grid_size),1)';
    index = index(1:length(channum));
    G = G(:,index);
    channum = channum(index);
end

if mod(theta,90) == 0
    index = reshape([1:prod(grid_size)]', grid_size(2), grid_size(1))';
    index = rot90(index,-theta/90)';
%     index = index(:);
%     index = reshape(index',prod(grid_size),1)';
%     index = index(1:length(channum));
%     G = G(:,index(index<=size(G,2)));
%     channum = channum(index(index<=size(G,2)));
    channum(index<=size(G,2)) = channum(index(index<=size(G,2)));
    plotme = index <= nchan;
    
%     G(:,max(index)) = 0;
    index(index>size(G,2) ) = [];
    G = G(:,index);
%     G(:,index>size(G,2)) = nan;
    nchan = max(index(:));
    if mod((theta/90),2) == 1
        grid_size = fliplr(grid_size);
    end
    theta = 0;
end
    


%if udflip == 1
%    I = [1:nchan]';
%    nusize = [ceil(sqrt(nchan)), ceil(nchan/sqrt(nchan))];
%    I = [I; zeros(prod(nusize) - nchan)];
%    I = reshape(I,nusize(2),nusize(1));
%    I = fliplr(I);
%    reshape(I,prod(nusize),1);
%    I = nonzeros(I)';
%    G = G(:,I);
%    channum = channum(I);
%end

%if theta ~=0
%    I = [1:nchan]';
%    nusize = [ceil(sqrt(nchan)), ceil(nchan/sqrt(nchan))];
%    I = [I; zeros(prod(nusize) - nchan)];
%    I = reshape(I,nusize(2),nusize(1));
%%    I = rot(I,theta/180*pi);
%    reshape(I,prod(nusize),1);
 %   I = nonzeros(I)';
%    G = G(:,I);
%    channum = channum(I);
%end

if isempty(chlbl)
    %chlbl = channum;
    prech = ' ';
    chlbl = [prech(ones(length(channum),1),:),num2str(channum')];
    chlbl = mat2cell(chlbl,ones(1,size(chlbl,1)),size(chlbl,2));
    %else
   % chlbl = chlbl(index);
end

axes_pos = [];
line_handles = [];

for i  = 1:nchan
   if ~plotme(i)
       continue
   end
    if new_plot
       subplot(grid_size(1), grid_size(2) ,i);
       axes_handles = [axes_handles;gca];
   else
       axes(axes_handles(i))
       hold on
    end

   if isempty(xdata)
       x = (1:size(G,1))./Fs - stonset;
   else
       x = xdata;
   end

   if i <= size(G,2)
       line_handles(end+1) = plot(x,G(:,i),plot_command{:});       
   end

   if ~isempty(linecolor)
       set(line_handles,'color',linecolor)
   end

   if length(axlim) == 4
       ax = axlim;
   elseif length(axlim) == 2
       ax(3:4) = axlim;
   end



   for j = 1:length(ax_command)
       eval(ax_command{j});
   end

   axis(ax)

   if flipy
       set(gca,'ydir','reverse')
   end


   if ~isempty(stonset)
       if length(stonset) > 1
         for sto = stonset
          line([sto sto],[ax(3) ax(4)],'linestyle',axis_line_style,'color',axis_line_color)
         end
       else
           line([0 0],[ax(3) ax(4)],'linestyle',axis_line_style,'color',axis_line_color)
%            Xdat = get(line_handles(end),'xdata');
%            set(line_handles(end),'xdata',Xdat - stonset);
%            currax = axis(gca);
%            axis(currax - [stonset, stonset, 0, 0]);
       end

   end

   if i<=length(channum) &&  ~ismember(channum(i),tickson) && ~strcmp(lower(tickson),'all')                                    % axis values showing for all channels in chlbls.
     set(gca,'XTickLabel','')
     set(gca,'YTickLabel','')
   elseif ~isempty(Fs)
       Xdat = get(line_handles(end),'xdata'); 
       set(line_handles(end),'xdata',Xdat*1000/Fs);
       currax = axis(gca);
       axis(currax.*[1000/Fs, 1000/Fs, 1, 1]);
       %L = get(gca,'xtick')*1000/Fs;
       %set(gca,'xticklabel',L);
%        xlabel('msec')
   end



   if ~isempty(horzline)
       currax = axis(gca);
       for hoz = horzline
           line([currax(1) currax(2)],[hoz hoz],'linestyle',axis_line_style,'color',axis_line_color)
       end
   end

   if ~iscell(chlbl)
       chlbl = num2cell(chlbl);
   end

   if length(chlbl) > 1 && i <= length(channum)
      %t = title(num2str(num2str(chlbl{i})),'fontsize', title_size);
       t = title(num2str(num2str(chlbl{channum(i)})),'fontsize', title_size);
   else
      t = title(num2str(num2str(chlbl{1})),'fontsize', title_size);
   end

    set(t,'units','normalized')
    tpos = get(t,'position');
    set(t,'position',[tpos(1),title_height,tpos(3)])
    %title(['Ch. ',num2str(chlbl(i))],'fontsize', 8)

   axes_pos = [axes_pos; get(gca,'position')];
end

if nargout > 0
    varargout(1) = {axes_handles};
end

if theta ~= 0           %rotates subplots around the center
   thrad = pi*theta/180;
   ROT = [cos(thrad) sin(thrad); -sin(thrad) cos(thrad)];                               %Rotation matrix

   cent = mean(axes_pos(:,1:2));                                %Center of rotation is center of axes locations. Change to [.5 .5] to make it the center of the figure.
   
   tot_width = max(axes_pos(:,1) + axes_pos(:,3)) -min(axes_pos(:,1));
   tot_height = max(axes_pos(:,2) + axes_pos(:,4)) -min(axes_pos(:,2));
  
   
   rot_adj = [abs(cos(thrad)), abs(sin(thrad))];                                            % adjustment so that axes lie within the original box
   pos_adjustment = [(rot_adj*[tot_width, tot_height]')/tot_height,  (rot_adj*[tot_height, tot_width]')/tot_width];  %
   
   axes_pos_cent = axes_pos(:,1:2) - cent(ones(1,length(axes_pos)),:);
   axes_pos_cent_adj = axes_pos_cent./pos_adjustment(ones(length(axes_pos_cent),1),:);

   axes_pos_rotated = (ROT*axes_pos_cent_adj(:,1:2)')' + cent(ones(length(axes_pos),1),:);
   nusize = axes_pos(:,3:4)./pos_adjustment(ones(length(axes_pos),1),:);
   new_axes = [axes_pos_rotated, nusize];
   
   for i = 1:length(axes_handles)
      set(axes_handles(i),'position',new_axes(i,:))
      axes_pos(i,:) = get(axes_handles(i),'position');    
   end
end
   
if backwhite   
    set(gcf,'Color',[1 1 1])
end

if exist('stim','var')
    axes;
    title(stim)
    axis off
end
