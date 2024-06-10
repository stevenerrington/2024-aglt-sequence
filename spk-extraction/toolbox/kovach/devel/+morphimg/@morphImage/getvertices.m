

function  getvertices(morphim,fixax,tesswgt)

%
%
% 



% X = morphim.X;
% Y = morphim.Y;
% VY = morphim.VY;
% VX = morphim.VX;
% tesswgt = morphim.tesswgt;

nim = length(morphim.images);

PtSelectionTol = 10;

fig = figure; 
for i = 1:nim
    colormap gray
    axvec(i) = subplot(1,nim,i); imagesc(morphim.images(i).X);
    axis image, axis off, hold on
end

morphim.tesswgt = tesswgt;
% axvec = [ax1 axvec(2)]; 

colors = {'b','g','r','c','m','y','w'};

%VX = [];
%VY = [];
p1 = [];
p2 = [];
j = 0;

quit = 0;
plotcolr = [];

plotopts= {'markersize',12};

pass2plot = {};
% pass2plot2 = {};

nv = max(morphim.images.nvert);
for i = 1:nim
   for k = 1:nv
%         
%         if k == currvertind
%             add = '+';
%         else
            add = '.';
%         end
%       
        pass2plot{1,i} =  morphim.images(i).vert(k,1);
        pass2plot{2,i} =  morphim.images(i).vert(k,2);
        pass2plot{3,i} =  [colors{mod(k,length(colors))+1},add];

        p1(end+1,i)  = plot(axvec(i),pass2plot{:,i},plotopts{:});
   end
end   
   
u = uicontrol('visible','on','units','normalized','position',[.46 .05 .1 .05],...
   'style','togglebutton','string','tesselate');
setappdata(fig,'buttonh',u);
% setappdata(fig,'axes',[ax1 axvec(2)]);
morphim.axes = axvec;
set(u,'Callback',@morphim.trplot)

% setappdata(fig,'VX',VX)
%  setappdata(fig,'VY',VY)
    
while ~quit
    
    colorind = mod(j,length(colors))+1;
    
    figure(fig)    
    refresh(fig)
    pause(.05)
    waitforbuttonpress
%     waitfor(fig,'userdata',0)
%     drawnow
    %     pause(1)
   
            
    %get(gcf,'selectiontype')
    if isempty(get(gcf,'currentcharacter'))
        quit = 0;
    elseif isspace(get(gcf,'currentcharacter'))
        break
    end
    
    currV = morphim.images(gca == axvec).vert;
    if isempty(currV)
        continue
    end

    currpt = get(gca,'currentpoint'); 
    currpt = currpt(1,1:2);
    
    if strcmp(get(gcf,'selectiontype'),'normal')
            
      try
        dv = ones(size(currV,1),1)*currpt - currV;
        dv = sqrt(sum(dv.^2,2));
        [mn,currvertind] = min(dv);
        if mn < PtSelectionTol
          if gca == axvec(1)
              pa = p1;
              pb = p2;
          elseif gca == axvec(2)
              pa = p2;
              pb = p1;
          end
          set(pa(currvertind),'xdata',currV(currvertind,1),'ydata',currV(currvertind,2),'marker','+')
          set(pb(currvertind),'marker','+')
          refresh(fig)
        else
          continue
        end    
      catch %#ok<*CTCH>
          continue
      end
    elseif strcmp(get(gcf,'selectiontype'),'alt')
       
       if fixax ~= 0, continue, end
       
 
       if fixax ~= 0, continue,end
 
       try 
        curraxdim = axis(gca);
        
        if ~((currpt(1) < 0) || (currpt(2) < 0) || (currpt(1) > curraxdim(2)) || (currpt(2) > curraxdim(4)))
         currV = [currV;round(currpt(1,:))]; %#ok<*AGROW>
         currvertind = size(currV,1);
         for i =1:nim
             morphim.images(i).vert = [morphim.images(i).vert;round(currpt(1,:))];
             p1(nv+1,i) = plot(axvec(i),morphim.images(i).vert(end,1),...
                 morphim.images(i).vert(end,2),[colors{colorind},'+'],plotopts{:});
         end
         
         plotcolr(end+1) = colorind;
         refresh(fig)
         j = j+1;
        else
         continue
        end
       catch
           continue
       end
    elseif strcmp(get(gcf,'selectiontype'),'extend')
       
       if fixax ~= 0, continue,end
       
       try
        dv = ones(size(currV,1),1)*currpt - currV;
        dv = sqrt(sum(dv.^2,2));
        [mn,currvertind] = min(dv);
       if mn < PtSelectionTol
          
%         currV(currvertind,:) = [];

        for i = 1:nim
            morphim.images(i).vert(currvertind,:) = [];
        end
        delete(p1(currvertind,:))
        p1(currvertind,:) = [];
        
         plotcolr(currvertind) = [];
        
        currvertind = [];
        
        refresh(fig) 
       else
           continue
       end
       catch caterr %#ok<NASGU>
           continue
       end
    else
%         currvertind = [];
	 	continue;   
    end
    
    if ~isempty(p1), set(p1,'zdata',1),end
%     if ~isempty(p2), set(p2,'zdata',1),end
    
   
%     morphim.VX = VX;
%     morphim.VY = VY;

    
    set(gcf,'selectiontype','normal')
    

    while 1
        
       isKey = waitforbuttonpress;
       %get(gcf,'selectiontype')
       
       if ~strcmp(get(gcf,'selectiontype'),'normal') || isKey
            break
       end
       
       if fixax ~= 0
         if axvec(fixax) == gca, continue, end
       end
       
       currpt = get(gca,'currentpoint');
       currpt = round(currpt(1,1:2));
       curraxdim = axis(gca);
       
       if ~((currpt(1) < 0) || (currpt(2) < 0) || (currpt(1) > curraxdim(2)) || (currpt(2) > curraxdim(4)))
            morphim.images(gca == axvec).vert(currvertind,:) = currpt;
            set(p1(currvertind,gca == axvec(1)),'xdata',...
                morphim.images(gca == axvec).vert(currvertind,1),...
                'ydata',morphim.images(gca == axvec).vert(currvertind,2))
       end

       
    end
    
     set(p1(currvertind,:),'marker','.')
%      
%     setappdata(fig,'VX',VX)
%     setappdata(fig,'VY',VY)
%     morphim.VX = VX;
% morphim.VY = VY;
% %     trh = getappdata(fig,'trh');
    trh=morphim.trh;
    if get(u,'value') 
        trplot(u)
    elseif all(ishandle(trh(:)))
        set(trh,'visible','off')
    end
end    

close(fig)