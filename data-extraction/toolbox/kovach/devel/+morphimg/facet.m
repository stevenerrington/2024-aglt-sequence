
classdef facet 
    
    %
    % class to store facet data
    %
    
    properties(SetAccess = private)
        
        trans = eye(3);
               
        idcode =0;
        
        vertex = [1 1 1]'*[0 0  1];
        reference = [1 1 1]'*[ 0 0 1];        
    end
    
    properties( Dependent = true)
        
        vert
        ref
    end
    
    methods 
        
        function me = facet(varargin)
                
             if nargin > 0
                  V = varargin{1};                      
                  
                  if nargin < 2 || isempty(varargin{2})
                      Vto = V;
                  else
                      Vto = varargin{2};
                  end
                  if nargin < 3
                      Tri = reshape(1:size(V,1),size(V,1)/3,3);
                  else
                      Tri = varargin{3};
                  end


                for i = 1:size(Tri,1)
                    me(i).vert = V(Tri(i,:),:);
                    me(i).ref = Vto(Tri(i,:),:);                
                end
                
            end
            
        end
        
        function me=set.vert(me,val)
            
            me.vertex(:,1:2) = val;
            me.trans = me.vertex^-1*me.reference;
        end

        function val = get.vert(me)            
            val = me.vertex(:,1:2);
        end
        
        
        function me = set.ref(me,val)
            
            me.reference(:,1:2) = val;
            me.trans = me.vertex^-1*me.reference;            
        end
        function val = get.ref(me)            
            val = me.reference(:,1:2);
        end
            
     
        function prd = mtimes(me,other)            
            
            
            for i = 1:length(me)
                prd(i) = me(i);
                prd(i).vert = other(i).vert;
                prd(i).ref = me(i).vert;
            end    
        end
        
        
        function d = det(me)  
            for i = 1:length(me)
                d(i) = det(me(i).trans);
            end
        end
        
        function isin = inside(me,xy)  
            isin = false(size(xy,1),length(me));
            for i = 1:length(me)
                vxy = [me(i).vert];
                
                isin(:,i) = inpolygon(xy(:,1),xy(:,2),vxy(:,1:2:end)',vxy(:,2:2:end)');
                
%                 me(i).ref = [0 0; 1 0; 0 1];
%                 xy(:,3) = 1;
%                 xytrans = xy*me(i).trans;
%                 isin(:,i) =  all(xytrans(:,1:2)>0,2) & sum(xytrans(:,1:2),2)<1;                
            end
        end
        
        function pp = patch(me,z)  
            xy = cat(2,me.vert);
            if nargin < 3 
                z = 1:length(me);
            end
            pp = patch(xy(:,1:2:end),xy(:,2:2:end),z); 
            
         end
            
    end
    
    
end