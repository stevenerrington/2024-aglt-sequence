
classdef helium < handle 
   
    properties (Dependent = true)
       status; 
    end
    properties
        
        job = '';
        description = ''
        
        email = '';
    end
    
    properties (SetAccess = private)
        
        path = '';
        compiled = false;
        ready = false;
        running = false;
        finished = false;
    end
    
    
    methods
        
        function me = helium
           
            a = uigetdir('C:\','Locate Helium mount point.');
            b = inputdlg('Which job (Enter to select existing)?');
            b = fullfile(a,'jobs',[b{:}]);
            if isempty(b)
                b = uigetdir(fullfile(a,'jobs','null','Choose job'));
            end
            [~,ft] = fileparts(fileparts(b));
            if ~strcmp(ft,'jobs')
                error('Must be a subdirectory of %s',fullfile(a,'jobs'));
            end
            if ~exist(b,'dir')
                mkdir(b)
                mkdir(fullfile(b,'in'))
                mkdir(fullfile(b,'out'))
                mkdir(fullfile(b,'bin'))
                mkdir(fullfile(b,'log'))
                mkdir(fullfile(b,'status'))
            end
            me.path = fullfile(b);
            
        end
        
        function compile(me)
            
            [fn,pth] = uigetfile('*.m','Select matlab script to compile');
            
            mccargs = {['-d ',fullfile(me.path,'bin')],'-m',fullfile(pth,fn)};
            mcc(mccargs{:}); % Matlab compiler
            me.compiled = true;
            me.logger(sprintf('#### Compiled binaries ####\nCommand: mcc %s\n#### )',sprintf('%s ',mccargs{:})));
        end
            
        function logger(me,msg,type)
            if nargin < 3
                type =0;
            end
            switch type
                case 0
                    logf = 'report.log';
                case 1
                    logf = 'status.log';
                case 2 
                    logf = 'error.log';
            end
            fid = fopen(fullfile(me.path,'log',logf),'a+');
            fprintf(fid,'\n%s: %s',datestr(now),msg);
            fclose(fid);

        end
        
        function notify(me,msg)            
            sendmail(me.email,sprintf('Notification for %s',me.job),msg);
        end
            
        
        function status
            
        end
    
    
    
        
    
    
    
end