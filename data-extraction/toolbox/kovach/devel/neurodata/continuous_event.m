classdef continuous_event < continuous
    
    %%% Class for continuous event data
    
    properties
              
      threshold = [];  
      minDelay = 0;
    end
    properties (Dependent = true)
        time
    end
    
  methods 
    function me = continuous_event(varargin)
        me = me@continuous(varargin{:});
    end
    
    function extract(me,side)
        if nargin<2
            side = 'start';
        end
        
        switch side
            case {'start','rise','rising'}
                trf = @(x)x;
            case {'end','fall','falling'}
                trf = @(x)-x;
            case {'both'}
                trf = @(x)abs(x);
        end
                    
        thr = me.data*sign(me.threshold) > abs(me.threshold);
        clwin = ones(round(me.minDelay*me.sampling_rate),1);
        clthr = convn(thr,clwin,'same')>0;
        clthr = cat(1,zeros(1,size(clthr,2)),clthr);
        evtind = diff(clthr)>0;
        
        
    end
  end
    

end