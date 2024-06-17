
classdef isolator_control < handle

%   iz2 = IZ2_controller   
%
% Creates a stimulator class object to control the IZ2 via a parallel port
% connection to the TDT system. This is designed to control stimulation 
% parameters and initiaon of the TDT IZ2 when the IZ2_control circuit is running. 
% 
% iz2.prime_stimulator : sends stimulation parameters to TDT and instructs
%                         the system to enable the stimulator.
% iz2.trigger          : Initialize stimulation.
%
% iz2 is an object with the following properties:
%
%     do_stimulate  : If set to 'off' or 0, iz2.trigger is disabled
%        stim_mode  : 'volt' or 1 - constant voltage
%                     'current' or 0 - constant current        
%            level  :  stimulation level in V for voltage or uA for
%                               current.
%      level_units : units of level 
%         duration : duration of pulse train in sec.
%        frequency : frequency of pulse train in Hz.
%         pos_chan : vector of positive polarity channels
%         neg_chan : vector of negative polarity channels
%         isprimed : ready to trigger if true
% 
%  
% 
%
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/TDT/@stimulator/stimulator.m $
% $Revision: 85 $
% $Date: 2012-03-03 13:57:49 -0600 (Sat, 03 Mar 2012) $
% $Author: ckovach $
% ------------------------------------------------



    properties ( Dependent=true )
        do_stimulate = 1;
        stim_mode= 0;
       %level = 0;
    end
     properties
    	
%         mode = 'ppt';
    %     mode = 'tdt';
%         level = 0;
        duration = 1;
        frequency = 50;
        pos_chan = 1;
        neg_chan = []; 
        level = 0;
     end
     
    
    methods
        
        function me = isolator_control(varargin)            
            me.pport = initializeLPT;
            initialize_stimulator(me,varargin{:});
        end
        
        function prime(me)
%             t1=tic;
            prime_stimulator(me);
%             toc(t1)
        end
        
        function varargout =trigger(me,code)
            if nargin < 2
                code = 1;
            end
             ppwrite = me.trigger_code*me.do_stimulate_val + me.no_stim_code*(1-me.do_stimulate_val);
            if me.isprimed
                me.LPT_pulse(ppwrite);

                if nargin > 1 && ~isempty(code)
                   me.putcode(code);
                end
                
            else
                error('Stimulator can''t be triggered until it is primed');
            end
            me.isprimed=false;
            if nargout > 0 
                varargout{1} = ppwrite;
            end
            
        end
        
        function varargout =putcode(me,code)
            %%% Put a code without altering stimulator settings
            %%% (verifies that the code is not 254)
            
%             if code > 128
%                 error('Code must be less than 128');
%             end
     
            if code ==254
                error('Code must not be 254');
            end

              
            me.LPT_pulse(code);
            if nargout > 0 
                varargout{1} = code;
            end
         
        end 
        
        function reset(me)
            me.LPT_pulse(me.reset_code);
            me.isprimed=0;
        end
        
      
        function set.stim_mode(a,b)
           switch lower(b)
               case {'0',0,'current','curr'}
                   a.stim_mode_val=0;
                   a.level_units = 'mA';
               case {'1',1,'voltage','volt'}
                   a.stim_mode_val = 1;
                   a.level_units = 'V';
               otherwise
                   error('Mode must be ''voltage'' or ''current''')
           end
                     
        end
        
        
        function b=get.stim_mode(a)
           switch a.stim_mode_val
               case 0
                   b='current';
               case 1
                   b='voltage';
           end
              
        end
        
        
         function b=get.do_stimulate(a)
           switch a.do_stimulate_val
               case 0
                   b='off';
               case 1
                   b='on';
           end
              
         end
         
         
         function set.do_stimulate(a,b)
           switch b
               case {'disabled',0,'off','no'}
                   a.do_stimulate_val=0;
               case {'enabled',1,'on','yes'}
                   a.do_stimulate_val=1;
           end
              
         end
       
          
         %%%
        function LPT_pulse(me,code)
            
            %%% Send a pulse over the lpt port of specified duration.
            me.pport.put(0)
            tic
            while toc < me.TTL_pulsewidth               
            end
            
            me.pport.put(code)
            tic
            while toc < me.TTL_pulsewidth               
            end
            me.pport.put(0)

          
          
        end
    end
    
    
    properties ( SetAccess = private )        
        isprimed = false;
        level_units ='uA';
    end
    properties ( Hidden = true, SetAccess = private )
        stim_mode_val = 0;
        do_stimulate_val = 1;
        level_val = 0;
        pport = 0;
        TTL_pulsewidth = .002; 
        reset_code = 255;
        trigger_code = 254;
        no_stim_code = 154;
    end
    properties ( SetAccess = private, GetAccess = private )        
        codeoffset =128;
   
       
%         wait = .01;
    end
end