classdef bismonitor < handle
    
    % Class for the interface with the coviden BIS monitor.
    %
    % bisobj = bismonitor;
    %
    % To record data and generate LPT evnt codes 
    %   
    %       bisobj.record([filename]);
    %
    % To stop
    %       
    %       bisobj.stop
    %
    %  LPT output encodes events as:
    %
    %     255--TIMESTAMP[0]--...--TIMESTAMP[7]--255--Ch12_BIS*2
    %
    %   where TIMESTAMP[k] is the kth of 8 bytes for the bismonitor
    %   timestamp encoded as datenum*1e6*24*3600 (matlab datestamp in
    %   microsecs)
    %
    properties
        serial_port = 'COM1';
        protocol = 'ASCII';
        recordFile = ''; 
%         lptfields = {'TIME','Ch3_BIS'}
%         get_data=false;
        poll_every = 1; % Polling rate every 1 sec. 
        header
        lpt_obj
        send_lpt_codes = true;
    end
    
    properties (Dependent = true)
       isConnected
       isRecording
    end
    properties (SetAccess = private)      
        serial_obj
        buffer =''
        fid
        numRecords=0;
        timer_obj = timer;
        bufferSize = 1024;
        parsedData
    end
    methods
       
        function me = bismonitor
            
            connect(me);
         
        end
        
        function connect(me)
            
            switch lower(me.protocol)
                case 'ascii'
                    BaudRate = 9600;
                    DataBits = 8;
                    StopBits = 1;
                    Parity = 'none'; 
                case {'binary','bin'}
                    BaudRate = 57600;
                    DataBits = 8;
                    StopBits = 1;
                    Parity = 'none';
                otherwise
                    error('Protocol must be ASCII or BINARY')
                               
            end
            if me.send_lpt_codes
                   me.lpt_obj = initializeLPT;
            end
            me.serial_obj = serial(me.serial_port,'BaudRate',BaudRate,'DataBits',DataBits,'StopBits',StopBits','Parity',Parity,'ReadAsyncMode','Continuous','inputBufferSize',me.bufferSize);
           fopen(me.serial_obj); 
           me.getHeader;
%            me.isConnected = strcmp(get(me.serial_obj,'Status'),'open')
        end
        
        function record(me,fname)
            % Start recording
            if nargin <2
                fname = me.recordFile;
            end
            while isempty(fname)
                fname = input('Please provide file name: ','s');

            end
            me.recordFile = fname;
            me.fid = fopen(fname,'a+');
            me.numRecords = 0;
%             me.isRecording = true;
%             me.get_data = true;
            me.flushBuffer;
            fprintf(me.fid,'\n--START RECORD--\n');
            out = query(me.serial_obj,'D');
            me.timer_obj = timer('Period',me.poll_every,'ExecutionMode','fixedRate',...
                        'TimerFcn',@(varargin)me.timerfcn,...
                        'StopFcn',@(varargin)me.stop);
             start(me.timer_obj)    
             if ~me.isRecording
                 warning('Failed to start recording!')
             end
%             while me.get_data
%                 me.buffer = fscanf(me.serial_obj);
%                 fprintf(me.fid,'%s',me.buffer);
%                 me.numRecords = me.numRecords+1;
%             end
%             me.isRecording = false;
%             fclose(me.fid);
        end
        function timerfcn(me,varargin)
                try
             fwrite(me.fid,deblank(me.readBuffer));
            if me.send_lpt_codes
                dat = me.parseBuffer;
                if isempty(dat)
                    return
                    warning('parseBuffer produced no output')
                end
                   %matlab time stamp in ms
                    timestamp = round(dat.TIME*3600*24*1e6);
                    timewords = bitand(round(timestamp.*2.^([0:-1:-8]*8)),255);

                    bisword = round(2*str2num(dat.Ch3_BIS));
                    me.lpt_obj.put(255);
                    pause = .01;
                    tic, while toc<=pause;end % much more accurate than pause
                    for k = 1:length(timewords)
                        me.lpt_obj.put(0);
                          tic, while toc<=pause;end
                        me.lpt_obj.put(timewords(k));
                          tic, while toc<=pause;end
                    end
                    tic, while toc<=pause;end
                    me.lpt_obj.put(0);
                  
                    me.lpt_obj.put(255);
                    tic, while toc<=pause*10;end    
                    me.lpt_obj.put(bisword);
                    tic, while toc<=pause;end
                    me.lpt_obj.put(0);
                    
                    fprintf(me.fid,'|TRUE');
             
            end
           
                catch err
                         warning('Skipped event with error %s',err.message)
                    fprintf(me.fid,'|FALSE');
                end
                  fprintf(me.fid,'\n');
        end
        function stop(me)
            fprintf(me.fid,'\n---STOP RECORD---');
            fclose(me.fid);
%             me.get_data = false;
        end
        
        function out = readBuffer(me)
            out ='';
            while me.serial_obj.BytesAvailable >0
                out = [out,char(fread(me.serial_obj,me.serial_obj.BytesAvailable))'];           
            end
            me.buffer = out;
        end
        
        function a = get.isConnected(me)
            a = strcmp(me.serial_obj.Status,'open');
        end
        function a = get.isRecording(me)
            a = strcmp(me.timer_obj.Running,'on');
        end
        function disconnect(me)
           fclose(me.serial_obj); 
%            me.isConnected = false;
        end
        
        function flushBuffer(me)
            n = 0;
            while me.serial_obj.BytesAvailable > 0
                fread(me.serial_obj,me.serial_obj.BytesAvailable);
                n=n+1;
                pause(.1)
            end
%             me.readBuffer;
            n=n+1;
            fprintf( '\n%i records flushed\n',n);
            
        end
        
        function out = parseBuffer(me)
            if isempty(me.header)
                me.getHeader;
            end
           out = [];
%             buff = me.readBuffer;
            buff = me.buffer;
            rxp = regexp(buff,'([^\n]*)','tokens');
            rxp = [rxp{:}];
            for k = 1:length(rxp)
                re = deblank(regexp(rxp{k},'[^|]*','match'));
                if isempty(re{1})
                    continue
                end
                try 
                    timestamp = datenum(re{1});
                    re{1} = timestamp;
                catch
                    continue
                end
                
%                  geti = ismember(me.header,me.parsefields);
                if isempty(re)
                    continue
                end
                re(end:36)=[];
                args = [me.header(1:length(re));re];
                out = struct(args{:});
                me.parsedData=out;
            end
                        
                
        end
        
        
        function out = getHeader(me)
          
           out = {'XX'};
           while ~strcmp(out{1},'TIME')
               me.flushBuffer;
                query(me.serial_obj,'C');
               out = query(me.serial_obj,'D');
               pause(1)
                   out = [out,me.readBuffer];
               out = regexp(out,'[^\n]*','match');
               out = regexp(out,'[^|]*','match'); 
               out = deblank(out{2});
           end
           [~,~,unqi] = unique(out,'stable');
           for k = 9:length(out)
               cs = cumsum(unqi==unqi(k));
               out{k}= strcat(sprintf('Ch%i_%s',cs(k),out{k}));
           end
           me.header = regexprep(out,'-','_');
        end
        
        
        function delete(me)
            disconnect(me); 
            delete(me.timer_obj);
            delete(me.serial_obj);
        end
    end
end