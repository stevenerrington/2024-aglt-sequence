classdef evnts < neurodata
    
    % Class for a series of discrete events
    
    properties (Dependent = true)
        
        time
        code
        message
        
    end
    
    properties 
        
        evnt = evnt([]);
        
    end
    
    methods
    
        function t = get.time(me)
            t = [me.evnt.time];
        end
        function a = get.code(me)
            a = [me.evnt.code];
        end
        function a = get.message(me)
            a = {me.evnt.message};
        end
        
        function set.time(me,t)
            for k = 1:length(t)
                me.evnt(k).time = t(k);
            end
        end
        function set.code(me,a)
            for k = 1:length(a)
                me.evnt(k).code = a(k);
            end
        end
        function set.message(me,a)
            for k = 1:length(a)
                me.evnt(k).message = a{k};
            end
        end
        
    end
end