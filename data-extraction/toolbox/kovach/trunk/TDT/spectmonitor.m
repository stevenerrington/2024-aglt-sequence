classdef spectmonitor < handle

    
% Defines a spectmonitor object. 
% Spectmonitor shows a time frequency sweep.
% Spectrogram is computed with a windowed mutlitaper.
%  
%
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/TDT/spectmonitor.m $
% $Revision: 128 $
% $Date: 2012-06-08 20:37:23 -0500 (Fri, 08 Jun 2012) $
% $Author: ckovach $
% ------------------------------------------------

    
  properties
 
      
        store= 'LFPx';
        
    
        tank = '';
        block = '';
        frqrange = [2 Inf];
    
        mS=0;
        S=[];
        ntaper = 2;
        track = false;
%         runtime=0;
        tstart=0;
     
   end
%%%%
    properties (SetAccess = private )
         RecIndx = 0;
           currtime = 0;
        tick = 0;
        currind = 0;   
        timerh = [];
          frq = [];
         tt=[]; 
        normalization = 'smoothed';
    end
%%%%%
   properties ( Hidden = true )
         W=[];
        V=[];
        t0=[];
        TTX=[];
        getfr = [];
        poverlap =[];
         Noverlap = 0;
        twinN=[];
        ftwinN=[];
          ttick =[];
        ttick0 = [];
         fig=[];
        im = [];
        pl = [];
        runstart = [];
        
        normalization_modes = {'smoothed','average','none'};
    
        ax = [];
        val_tstep = .05;
        val_twin = 10;
        val_fftwin = 1;
        val_fs = 1e3;
        use_listdlg = true;
        val_speed = 1;
        val_runtime = 0;
   end
   
%%%%
    properties (Dependent=true)
        running
        channel
        fs 
        tstep
        twin
        fftwin
        speed
        runtime
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%      
    methods
        
        function me = spectmonitor(varargin)            
            me.TTX=actxcontrol('TTank.x');
            invoke(me.TTX,'ConnectServer','Local','Me');
%            
                me.gettank;
                me.getblock;
%             me.TTX.SetGlobalV('Channel',me.channel);
            stores = {'LFPx','PDes','Inpt'};
            st = listdlg('ListString',stores);
            
            me.store = stores{st};
            
            me.TTX.SetGlobalStringV('Options','FILTERED');
           
            me.timerh = timer('Period',me.tstep,'ExecutionMode','fixedRate',...
                'TimerFcn',@(a,b)me.update_figure);
            
             me.poverlap = (me.fftwin-me.tstep)./me.fftwin;

         
            me.fs = 1e3;
            
            me.initializeS;
           
            me.fig = figure;
            me.im = imagesc(me.tt,me.frq,me.S);
            me.ax = gca;
            hold on;
            me.pl = plot(0,0);

            set(me.fig,'DeleteFcn',@(a,b,c) delete(me));

            me.channel=0;
            
%             [me.W,me.V] = dpss(me.ftwinN,me.ntaper/2);
           
            axis([0 me.twin me.frq([1 end])]);
        end
%%%%
        function set.channel(me,channel)
            if channel==0
                channel = inputdlg('Which channel? ');
                channel = str2double(channel{1});
            end
            me.TTX.SetGlobalV('Channel',channel);
            title(me.ax,sprintf('Channel %i',channel));
        end
%%%%
        function ch= get.channel(me)
            ch=me.TTX.GetGlobalV('Channel');
        end
%%%%
        function set.fs(me,fs)
            
            me.TTX.SetGlobalV('WaveSF',fs);
            me.val_fs = fs; 
            me.initializeS;
        end
%%%%        
        function x = get.fs(me)
            %%% Set Sampling freq
            x = me.TTX.GetGlobalV('WaveSF');
        end
%%%%
        function x=get.tstep(me)
            x = me.val_tstep;
        end
%%%%
        function x=get.twin(me)
            x = me.val_twin;
        end
%%%%
        function x=get.fftwin(me)
            x = me.val_fftwin;
            
        end
%%%%
        function x=get.speed(me)
            x = me.val_speed;            
        end        
%%%%
        function  set.speed(me,speed)
            me.val_speed=speed;
            isr = isequal(me.timerh.Running,'on');
            if isr, stop(me.timerh); end
            set(me.timerh,'Period',me.tstep./speed);
            if isr, start(me.timerh); end
        end
        
%%%%
        function set.tstep(me,x)
            me.val_tstep=x;
            me.initializeS;
        end
%%%%
        function set.twin(me,x)
            me.val_twin=x;
            me.initializeS;
        end
%%%%%
        function set.fftwin(me,x)
            me.val_fftwin=x;
            me.initializeS;
        end
%%%%%
    function set.runtime(me,x)
            me.val_runtime=x;
            me.runstart = now*24*3600;
            me.tstart= x;      
            me.track = 0;
            me.tick = 0;
    end
%%%%
          function x=get.runtime(me)
            x=me.val_runtime;
        end
%%%%%
        function normalize(me)           
            %%% Choose normalization type
            inp =  listdlg('Liststring',me.normalization_modes,'selectionmode','single');            
            me.normalization = me.normalization_modes{inp};
            
        end
 %%%%%
         function run(me)
             
             %%% Start plotting
           me.TTX.CreateEpocIndexing();
           if me.channel == 0
               warning(sprintf('Channel set to 0! Are you sure you want to continue?\nCtrl+C to stop.'))                %#ok<WNTAG,SPWRN>
           end
           me.runstart=now*24*3600;
            start(me.timerh);
            me.tick=0;
         end
 %%%%%
         function stop(me)
             
             	me.TTX.CreateEpocIndexing();
                stop(me.timerh);
             
         end
%%%%
         function running = get.running(me)
            %%% Show if running
             running = me.timerh.Running;
         end   
%%%%
        function delete(me)
            me.TTX.ClearEpocIndexing;
            me.TTX.CloseTank;
            me.TTX.ReleaseServer;  
            
            delete(me.timerh)
            set(me.fig,'DeleteFcn',@(a,b,c)fprintf('\nGoodbye...\n'));

            close(me.fig)
        end
%%%%%%     
        function update_figure(me)
            
            
           %%% Update the figure 


            me.currind = mod(me.tick,me.twinN)+1;
            me.tick = me.tick+1;
            me.TTX.ResetFilters;  
             
%             me.runtime = me.speed*toc(me.runstart);
             me.val_runtime = me.speed*(now*24*3600-me.runstart)+me.tstart;
            if me.track
                 me.RecIndx=me.TTX.ReadEventsV(1e6,me.store,me.channel,0,0,0,'NEW')+me.RecIndx; 
            else
                 me.RecIndx=me.TTX.ReadEventsV(1e6,me.store,me.channel,0,me.runtime-me.fftwin/2,me.runtime+me.fftwin/2,'ALL'); 
            end
            ts = invoke(me.TTX, 'ParseEvInfoV', 0, 100, 6);
            ts = ts(ts>0);
            if me.currind == 1;
                me.t0=round(ts(end)*10)/10;
            end

%             geti=find(diff(ts(end)-ts <= me.fftwin)==1);

            me.TTX.SetGlobalV('T1',ts(end)-me.fftwin);
            me.TTX.SetGlobalV('T2',ts(end));
        %     x =  me.TTX.ReadWavesOnTimeRangeV(me.store,me.channel);
            x =  me.TTX.ReadWavesV(me.store);

            fx = log(abs(fft(repmat(x(:,end),1,me.ntaper).*me.W)).^2*me.V);
        %     [S,sf,st] = specgram(x,1048,d.fs(1));

            me.ttick(me.ttick0<me.tt(me.currind))=me.ttick0(me.ttick0<me.tt(me.currind))+me.t0;



            me.mS = (me.mS*me.tick+fx(me.getfr))/(me.tick+1);   
            switch me.normalization
                case 'none'
                    s = fx(me.getfr);
                case 'smoothed'
                    b = polyfit(log(me.frq(me.frq>0))',me.mS(me.frq>0),6);
                    pv = polyval(b,log(me.frq)');
                    s = fx(me.getfr) - pv;
                case 'average'
                    s = fx(me.getfr)-me.mS;
                    
            end

            me.S(:,me.currind)= s;
            set(me.im,'CData',me.S,'ydata',me.frq([1 end]))
            set(me.ax,'xticklabel',me.ttick);
            set(me.pl,'xdata',me.tt(me.currind)*[1 1] + .5*diff(me.tt(1:2)),'ydata',me.frq([1 end]))
            drawnow

%             t=now;
        end
%%%%%%
        function initializeS(me)
            %%% Create S of the appropriate size and shape
            
             me.twinN = ceil(me.twin./me.tstep);
             me.ftwinN = me.fftwin*me.fs;

            me.Noverlap = round(me.poverlap*me.ftwinN);
            
            freq = (0:me.ftwinN-1)./me.ftwinN*me.fs;
            me.getfr = freq <= me.fs/2 & me.frqrange(1)<freq & me.frqrange(2)>freq;
            me.frq = freq(me.getfr);
            me.tt = (0:me.ftwinN - me.Noverlap:me.twin*me.fs)./me.fs;
 
            me.ttick = linspace(me.tt(1),me.tt(end),11);
            me.ttick0 = me.ttick;
        
            me.S = nan(sum(me.getfr),me.twinN);
            
            [me.W,me.V] = dpss(me.ftwinN,me.ntaper/2);
            
            me.mS=0;
            
        end
        
    
        function gettank(me,tankname)

%             use_listdlg = true;

            if  nargin < 2 || isempty(tankname)
               inp = 0; 
               tnames = me.tdtTanks;
               ntanks = length(tnames);
               if nargin == 1|| isempty(tankname) 
                   
                   [~,sind] = sort(lower(tnames));
                   tnames = tnames(sind);
                   tnames{end+1} = 'Other (Not Listed)';

                   if me.use_listdlg && usejava('jvm')
                       inp = listdlg('PromptString','Which Tank?',...
                                     'ListString',tnames,'selectionmode','single');
                   else

                       inpstr = sprintf('%i. %%s\n',1:length(tnames));
                       inp = input(sprintf(['Which Tank?\n',inpstr,':'],tnames{:}));
                       fprintf('\n\n\n\n\n')
                   end
                   if isempty(inp)
                       return
                   end
                   tankname = tnames{inp};
               end
               if ~ismember(tankname,tnames(1:ntanks)) || inp > ntanks 
                fprintf('\nOpening unregistered tank. Please select data tank file.\n\n\n')
                [tankname,tpath] = uigetfile({'*.tev','TDT Data Tank'});
                tankname = strtok(tankname,'.');

               end

              if ~ismember(tankname,tnames(1:ntanks))
               invoke(me.TTX,'AddTank',tankname,tpath);
              end

            end
             me.tank = tankname;
             me.TTX.OpenTank(me.tank,'r');
        end
%%%%
         function getblock(me,blocknames)

%             use_listdlg = true;
            i = 0;
            if nargin < 3 || isempty(blocknames)
                while 1
                   i = i+1; 
                    blks{i} = invoke(me.TTX,'QueryBlockname',i); %#ok<*AGROW>
                    if isempty(blks{i})
                        blks(i) = [];
                        break
                    end

                end

                if length(blks) > 1
                    if me.use_listdlg && usejava('jvm')
                        inp =listdlg('PromptString','Which block?',...
                                    'ListString',blks,'selectionmode','single');
                    else
                        inpstr = sprintf('%i. %%s\n',1:length(blks));
                        inp = input(sprintf(['Which Block?\n',inpstr,':'],blks{:}));
                    end
                    if isempty(inp)
                        return
                    end
                    blocknames = blks(inp);
                else
                    blocknames = blks(1);
                end
            else
                 if isfloat(blocknames)
                    for q = 1:length(blocknames)
                       blockn{q} = invoke(me.TTX,'QueryBlockName',blocknames(q));
                       if isempty(blockn)
                           error('There aren''t that many blocks.')
                       end
                    end
                    blocknames = blockn;
                end


            end
            me.block = blocknames{1};
            me.TTX.SelectBlock(me.block);
            
            if isequal(lower(me.block),'tempblk')
                me.track = true;
            end
            
        end
        
%%%%
        function T = tdtTanks(me)
            %
            %Generates a list of registered tdt tanks

            %me.TTX=actxcontrol('TTank.x');
            %invoke(me.TTX,'ConnectServer','Local','Me');

            i = 0;
            T = {''};
            while 1

                tname = invoke(me.TTX,'GetEnumTank',i);
                if isempty(tname)
                    break
                end
                T{i+1} = tname;
                i = i+1;
            end
        end

    end
end