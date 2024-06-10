classdef evnts < neurodata
    
    
    
    properties (Dependent = true)
        
        time
        code
        message
        
    end
    
    properties 
        
        evnt
        
    end
    
    methods
    
        function t = get.time(me)
            t = [me.evnt.evtime];
        end
        function a = get.code(me)
            a = [me.evnt.evcode];
        end
        function a = get.message(me)
            a = {me.evnt.message};
        end
        
        function set.time(me,t)
            for k = 1:length(t)
                me(k).evnt(k).time = t(k);
            end
        end
        function set.code(me,a)
            for k = 1:length(a)
                me(k).evnt(k).code = a(k);
            end
        end
        function set.message(me,a)
            for k = 1:length(a)
                me(k).evnt(k).message = a{k};
            end
        end
        
    end
    end