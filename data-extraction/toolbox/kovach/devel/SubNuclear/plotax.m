
classdef plotax < handle
% Class definition of the plotax object for volume view. It contains
% information about the plane of view for a given axis.

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/plotax.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

% C. Kovach 2013
    properties
        label = '';
        notes= '';
      
        h = -1;
        objectid = -1;
        handles;

%         normvec;
%         axvec;
        rot = 0;
        transpose = 0;
        parent;
        axdim;
        sisters;
        
    end
     properties (Hidden = true)
         basetr = eye(3,2);
         normv;
     end
     properties (Dependent = true)
         basetrmat;
         normvec;
     end
     
    methods
        
        function me = plotax(varargin)
            
            if nargin > 0 &&  isa(varargin{1},'plotax')
                
                props = properties(varargin{1});
                for i = 1:length(props)
                   me.(props{i}) = varargin{1}.(props{i}); 
                end
                
            elseif nargin > 0

                me.h = varargin{1};
            end
            
        end
        %%%
        
        function set.basetrmat(me,a)
            me.basetr = a;
            for i = 1:length(me.sisters)
                me.sisters(i).basetr= a;
            end
        end
        %%%
        function b=get.basetrmat(me)
            b = me.basetr;
        end
            %%%
        
        function set.normvec(me,a)
            me.normv = a;
            for i = 1:length(me.sisters)
                me.sisters(i).normv= a;
            end
        end
        %%%
        function b=get.normvec(me)
            b = me.normv;
        end 
        %%%
        function ish = isopen(me)
            ish = ishandle([me.h]);
        end
        %%%
        function c = children(me)
           c = get(me.h,'children'); 
        end
        function out= get(me,varargin)
            out = get(me.h,varargin{:});
        end
        function set(me,varargin)
            set(me.h,varargin{:});
        end
    
    %%%
         function ptax  = vol2ax(ax,ptvol)
       % Project volume coordinate into axis space
          ptax = [ptvol, ones(size(ptvol,1),1)]*ax.trmat;
     end
      %%%
     function ptvol = ax2vol(ax,ptax)
        
         % Project axis coordinate into volume space
        vv = ax.parent;
        crp = [vv.current_point 1];
        vc = size(vv.volumes(vv.current.volumes==[vv.volumes]).image.Data)/2;
       
        T = ax.trmat;
       nv = [ax.normvec,0];
       A = [T,nv'];
         A(4,4) = 1;
%         axcp = get(ax,'CurrentPoint');
%           
%           if isempty(ax.pivot)
%               ax.pivot = (crp(1:3)-vc)*nv(1:3)'*nv(1:3)+vc;
%           end
          
          p0 = [(crp(1:3)-vc)*nv(1:3)'*nv(1:3)+vc,1];
        newcp = [ptax - repmat(p0,size(ptax,1),1)*T,zeros(size(ptax,1),2)]*A^-1+repmat(p0,size(ptax,1),1);
%         crpnew = crp - (crp*axv')*axv + axcp(1,[2 1])*axv;
        ptvol = newcp(:,1:3);
         
     end
     function trbutton(ax,src)
            
            
            switch src
                case ax.handles.transpose % transpose
                   ax.transpose = get(src,'value'); 
                case ax.handles.rotate % Rotate 90 degrees
                    ax.rot=mod(ax.rot+1,4);
            end
              ax.parent.resetaxis(ax);
           
     end
     function T=trmat(ax)
            %Return the rotated transformation matrix 

            axd = ax.axdim;
            axd(5) =1;

           Ttranspose = [0 1 0 ; 1 0 0; 0 0 1];
           Trot = eye(3);
           if ax.rot>0
               
               Trot1 =  axd([1 3 5; 1 4 5; 2 4 5; 2 3 5]) \axd([ 1 2 5; 4 2 5; 4 3 5; 1 3 5]);
                   Trot = Trot1;
               if ax.rot > 1
                   
                   Trot2 =  axd([ 1 2 5; 4 2 5; 4 3 5; 1 3 5]) \axd([ 2 4 5; 2 3 5;1 3 5; 1 4 5]);
                   Trot = Trot*Trot2;
                   
                   if ax.rot >2
                       Trot = Trot*Trot1;

                   end
               end
           end
               
           tm=ax.basetrmat;
           tm(4,3) = 1;
           T = tm*(Ttranspose^ax.transpose)*Trot;
           T = T(:,1:2);
           
     end
 end
end