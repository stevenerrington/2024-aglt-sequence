classdef continuous_data < continuous
    
    %%% Class for continuous neuro data
    
    properties
        number
        contact
    
    end
  
    methods 
        function me = continuous_data(varargin)
            me = me@continuous(varargin{:});
        end
    end
end