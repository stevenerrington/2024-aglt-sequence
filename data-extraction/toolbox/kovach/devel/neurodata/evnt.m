classdef evnt 
    
        
    % Class for a single event

    properties
        
        time = nan;
        code = nan;
        message = '';        
    end
    methods
        function me = evnt(varargin)
            if nargin>0 && isempty(varargin{1})
                me = me([]);
            elseif nargin >0
                me = repmat(me,varargin{:});
            end
        end
    end
end