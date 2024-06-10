% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/transforms.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

classdef transforms   
    %Class definition for a volumeview transform class.
    properties
        label = [];
        notes= [];
        file = '';
        units = 'units';
        type = 'linear';
        objectid = -1;
        tr
        itr
    end
    properties (SetAccess=private)
        mytr = [];
        chain = [];
    end
    properties (Dependent=true)
        trmat = [];        
    end
        
   methods
        function me = transforms(varargin)
            
             for i = 1:2:length(varargin)
                me.(varargin{i})=varargin{i+1};
            end
        end
        
        function  a = invtrmat(me)
            if ~isempty(me.trmat)
                a = me.trmat^-1;
            else
                a = [];
            end
        end
        
        function me = set.trmat(me,T)
            if size(T,2)< size(T,1)
                T(end,end+1) = 1;
            end
            me.mytr = T;
            me.chain = me;
            me.tr = @(pt)me.lintr(pt);
            me.itr = @(pt)me.linitr(pt);
       
        end
       function T = get.trmat(me)
             T=me.mytr ;
        end

     
        function b = mtimes(me,a)
            b = transforms('label',sprintf('%s*%s',me.label,a.label),'trmat',me.trmat*a.trmat,'chain',[me.chain,a.chain]);
        end
        function b = eq(me,a)
           if length(me)==1
               x1 = me.trmat;
               c = a;               
           elseif length(a)==1
              x1 = a.trmat;
               c = me;
           else
               b = arrayfun(@(d,e) numel(d.trmat)==numel(e.trmat) && all(d.trmat(:)==e.trmat(:)),me,a);
               return
           end
           b = arrayfun(@(d) numel(d.trmat)==numel(x1) && all(d.trmat(:)==x1(:)),c);

        end
        function me = compute(me,fromx,tox)
            
            nd = size(fromx,2);
            switch nd
                case 3
                    fromx(:,4) = 1;
                case 2
                     tox(:,3) = 1;
            end
            me.trmat = tox\fromx;
            
        end
 
       function a = lintr(me,pt)
             
                % apply transform
                nd = size(pt,2);

                if nd < size(me.trmat,1)
                    pt(:,end+1) = 1;
                end
                a = pt*me.mytr;
                a = a(:,1:end-1);
       end
       function a = linitr(me,pt)
         % apply inverse transform
                nd = size(pt,2);
                if nd < size(me.trmat,1)
                    pt(:,end+1) = 1;
                end
                a = pt*me.mytr^-1;
                a = a(:,1:end-1);
       end
   end
end