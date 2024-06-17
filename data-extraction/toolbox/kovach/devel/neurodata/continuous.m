classdef continuous < neurodata
    
    %%% Class for continuous data.
  properties
      
        channel
        sampling_rate
        data
    
    end
    properties (Dependent = true)
        time
        nsamples
        nchannels    
    end
    
     methods
       
        function me=continuous(varargin)
            if nargin > 0 && isnumeric(varargin{1})
                me.data = varargin{1};
                varargin(1) = [];
            end
            i = 1;
            fn = fieldnames(me);
            while i < length(varargin)
                switch lower(varargin{i})
                    case 'fs'
                       me.sampling_rate =varargin{i+1};
                    otherwise
                    if ~ismember(varargin{i},fn)
                        error('%s is not a valid property.',varargin{i})
                    end
                    me.(varargin{i}) = varargin{i+1};
                end
                i = i+2;
            end
        end
    
         function tt = get.time(me)
           tt= (0:size(me.data,1)-1)./me.sampling_rate;
        end
        function ns = get.nsamples(me)
           ns = size(me.data,1);
        end
        function nch = get.nchannels(me)
          nch= size(me.data,2);
        end
        
        function sr = fs(me)
            sr = me.sampling_rate;
        end
        
     end
end