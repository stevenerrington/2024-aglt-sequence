
classdef imagedat
    
    %
    % class to store image data
    %
    
    properties 
        
        X = [];
        
        facets;
        vert = [];
        mapto=[];
        notes ='';      
       
        detA = [];
    end
    
    methods 
        
        function me = imagedat(X,VX)
            
            if nargin < 2
                VX = [];
            end
            if nargin < 1
                X = [];
            end
            me.X = X;
            me.vert = VX;
            
        end
        

        function sz = imsize(me) 
            for i = 1:length(me)
                sz(i,:) = [size(me(i).X,1),size(me(i).X,2)];
            end
        end
        
        function tc=truecolor(me)
            for i = 1:length(me)
                tc= size(me(i).X,3)==3;
            end
        end
        
        function f = fmt(me)
            f = class(me.X);  
        end
        
        function rg=range(me)
            for i = 1:length(me)
                rg(i,:) = [min(me(i).X(:)) max(me(i).X(:))];
            end
        end
        
        function x=rgb(me,lev)
           
            if ~me.truecolor
                x=cat(3,me.X);
            else

                switch lev
                    case {1,'r','R'} 
                        indx = 1;
                    case {2,'g','G'}
                        indx = 2;
                    case {3,'b','B'}
                        indx = 3;
                    otherwise 
                        indx = [];
                end
                for i = 1:length(me)
                    x(i,:) = me(i).X(:,:,indx);
                end
            end
        end
        
        function nv = nvert(me)
            for i = 1:length(me)
                nv(i) = size(me(i).vert,1);
            end
        end
        function X = double(me)
            
            if ~ismember(me.fmt,{'double','single',''})
                imx = double(intmax(me.fmt));
                X = double(me.X)./imx;
            else
                X= double(me.X);
            end
            
        end    
        
        function X = uint8(me)
            
            if ismember(me.fmt,{'double','single'})
                imx = double(intmax('uint8'));
                X = uint8(me.X*imx);
            else
                X=uint8(me.X);
            end
            
        end    
        
        
        %%% Weighted sum over the vertices
        function Vtot=weightedsum(me,wgt)
            if nargin < 2
                wgt = ones(size(me))./length(me);
            end
            Vtot=0;
            for i = 1:length(me)
                Vtot = wgt(i)*me(i).vert+Vtot;
            end
        end
            
        
    end
end