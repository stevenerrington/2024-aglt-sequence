
classdef fitgrid < SubModules.module

%
% Fitgrid interpolates a grid of specified size using the selected points 
% as control points. 
%
% It accomplishes this with either linear and non-linear warping.
%
% The default method is linear if 4 or fewer points are given (a minum of 3
% is reuired) and non-linear thin-plate spline warping otherwise.
%
% see also FITGRID


% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/+SubModules/sliceview.m $
% $Revision: 287 $
% $Date: 2013-07-29 14:39:16 -0500 (Mon, 29 Jul 2013) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 2013
    properties
        usedefault = true;
        default_size = [8 12];
        default_order = [3 3; 1 3; 3 2; 1 2]; 
        transf;
        points;
    end
    properties (SetAccess = private, Hidden = true)
        fig;
        szh;
        tabh;  
        buth;
        warph; 
        lblh;
    end
    properties (Dependent = true)
        gridpoints
        size    
        ptlabel
    end
    methods


        function me = fitgrid(varargin)   

            me.initialize(varargin{:});
            vv = me.parent;
            pts = vv.current.points;
           
            me.points=pts;

%             lbls = {pts.label};

            
            me.fig = figure('Name','Grid interpolation','numbertitle','off');
            set(me.fig,'units','characters','position',[60   30   50   30])
            set(me.fig,'DeleteFcn',@(a,b)delete(me))
           
        %                 cpt(1) = uicontrol('style','list','max',2,'units','normalized','position',[.1, .1, .1,.8],'string',lbls);
             sz = uitable('columneditable',true,'units','normalized','rowname','grid size','columnname',{'I','J'});
               pos = [.1, .8 0 0] + get(sz,'extent')*[1 0 0 0; 0 0 0 -1; 0 0 1 0; 0 0 0 1]';
            set(sz,'position',pos,'celleditcallback',@(src,evnt)me.szcallback(src,evnt),'data',me.default_size)
             me.szh = sz;

              tab = uitable('columneditable',[false,true(1,2)],'units','normalized','data',num2cell(zeros(length(pts),3)),'rowname',{});
              set(tab,'cellselectioncallback',@(src,evnt)me.tabselection(src,evnt))
              pos = [.1, .6 0 0] + get(tab,'extent')*[1 0 0 0; 0 0 0 -1; 0 0 1 0; 0 0 0 1]';
              pos(2) = max(pos(2),.1);
             set(tab,'position',pos,'celleditcallback',@(src,evnt)me.tabcallback(src,evnt),'columnname',{'Point','Grid I','Grid J'})             
             me.tabh = tab;
              
             me.buth = uicontrol('style','pushbutton','units','normalized','position',[.6, .85, .2,.1],'string','OK','callback',@(src,evnt)okcallback(me,src,evnt));
             me.warph(1) = uicontrol('style','radio','units','normalized','position',[.3, .875, .2,.06],'string','Linear','callback',@(src,evnt)wrpcallback(me,src,evnt));
             me.warph(2) = uicontrol('style','radio','units','normalized','position',[.3, .825, .2,.06],'string','TPS','callback',@(src,evnt)wrpcallback(me,src,evnt));
             uicontrol('style','text','units','normalized','position',[.1, .875, .2,.05],'string','Label','callback',@(src,evnt)wrpcallback(me,src,evnt));
              me.lblh =uicontrol('style','edit','units','normalized','background',[1 1 1],'position',[.1, .825, .2,.05],'string','Grid %i');
            
             if length(pts)<=4
                set(me.warph(1),'value',1);
                me.wrpcallback(me.warph(1)) 
             else
                set(me.warph(2),'value',1);
                me.wrpcallback(me.warph(2)) 
             end            
             me.szcallback();


        end
        function tabselection(me,~,evnt)
           if evnt.Indices(2) == 1
              me.parent.current_point = me.points(evnt.Indices(1)).coord;  
           end
        end
        function a = get.size(me)
           a =get(me.szh,'data'); 
        end
        function a = get.gridpoints(me)
           X =get(me.tabh,'data'); 
           a = cell2mat(X(:,2:end));
        end
         function a = get.ptlabel(me)
           a =get(me.lblh,'string'); 
        end
        function update(me) %#ok<MANU>
            % Nothing to update
        end
        %%%
        function szcallback(me,~,~)
            if me.usedefault
                X = get(me.szh,'data');
                X(3) = 1;
                X(4) = 0;
                txmat = me.default_order;
                txmat = txmat(1:min(length(me.points),4),:);
                txmat(min(length(me.points),4)+1:length(me.points),:) = 4;
         
                 set(me.tabh,'data',cat(2,{me.points.label}',num2cell(X(txmat))));
            end
        end
        
        function tabcallback(me,~,~)
           me.usedefault = false; 
        end
        function wrpcallback(me,src,~) 
            val = get(me.warph(me.warph==src),'value');
            set(me.warph(me.warph~=src),'value',1-val);
        end
        
        function okcallback(me,~,~) 
                       
           pts = me.points;
           Y = cat(1,pts.coord);
%            X = get(me.tabh,'data');
           warptype = find(cell2mat(get(me.warph,'value')));
         
           switch warptype
               case {1,2}
                 pts = fitgrid(Y,me.size,me.gridpoints,warptype==1,me.ptlabel);     
           end
           me.parent.addpoint(pts);
           if ~strcmp(get(me.fig,'beingdeleted'),'on')
               delete(me.fig);
           end
        end
    end
end