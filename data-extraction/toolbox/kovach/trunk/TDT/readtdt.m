 function varargout=readtdt(tankname,blocknames,DataCode,chnum,TT,FS,dryrun,interactive)


% D =readtdt(tankname,blockname,DataCode,chnum)
%
% Returns  structure D with fields:
% 
% D.dat - data
% D.fs  - sampling frequency for each sweep
% D.t   - onset time in seconds for each sweep
% D.chan - channel number for each sweep
% D.tankname - name of source tank
% D.block  - name of source block
% D.code   - name of source store
%
% If any or all arguments are empty, a menu of possible choices 
% for the corresponding argument will appear.
%
% [D,ttank] =readtdt(tankname,blockname,DataCode,chnum)
%
% Returns both D and the open activeX control ttank, which is otherwise destroyed
% at the termination of the script.
%
%
% D =readtdt(tankname,blockname,DataCode,chnum, ttank)
%
% Uses the aleady open control ttank rather than creating a new one.
%
% Reusing ttank (rather than allowing a new control to be created and destroyed
% every time readtdt is called) makes things faster, reduces memory
% leaks and is highly advisable if calling readtdt multiple times.
% 
% D =readtdt(tankname,blockname,DataCode,chnum, ttank,FS)
%
%  FS - specifies alternate sampling frequency at which to withdraw data.
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/TDT/readtdt.m $
% $Revision: 978 $
% $Date: 2018-02-08 16:55:41 -0600 (Thu, 08 Feb 2018) $
% $Author: ckovach $
% ------------------------------------------------

%Written by C Kovach, 2006
%
% 5/29/07:  modified to return any event store, not just Evnt
%           and to return all channels as columns when chnum = 0.
% 9/24/10:  specify sampling frequency at which to extract data.
% 6/26/12:  A bug in tdt extraction software causes epoc codes to be missed
%           with older data tanks. Now checking for this error and
%           extracting data through parseEvV


use_listdlg = true; %%% Use GUI interface if true and java is enabled


if nargin < 5 || isempty(TT)
    TT=actxcontrol('TTank.x');
    invoke(TT,'ConnectServer','Local','Me');
end
if nargin < 8 || isempty(interactive)
    interactive = true;
end
tpath='';    
if  nargin < 1 || isempty(tankname)
   inp = 0; 
   tnames = tdtTanks(TT);
   ntanks = length(tnames);
   if (nargin == 0 || isempty(tankname))&& interactive  
       
       [~,sind] = sort(lower(tnames));
       tnames = tnames(sind);
       tnames{end+1} = 'Other (Not Listed)';
          if use_listdlg && usejava('jvm')
                       inp = listdlg('PromptString','Which Tank?',...
                                     'ListString',tnames,'selectionmode','single');
           else

               inpstr = sprintf('%i. %%s\n',1:length(tnames));
               inp = input(sprintf(['Which Tank?\n',inpstr,':'],tnames{:}));
               fprintf('\n\n\n\n\n')
           end
       tankname = tnames{inp};
   else
       tankname = '';
   end
   if ~ismember(tankname,tnames(1:ntanks)) || inp > ntanks 
    fprintf('\nOpening unregistered tank. Please select data tank file.\n\n\n')
    [tpath,tankname,ext] = fileparts(tankname);
    if isempty(tpath)
        [tankname,tpath] = uigetfile({'*.tev','TDT Data Tank'});
    end
    tankname = strtok(tankname,'.');
     
   end
  
  if ~ismember(tankname,tnames(1:ntanks))
   invoke(TT,'AddTank',tankname,tpath);
  end

end

if nargin <= 1 || isempty(blocknames)
    blocknames = {};
end


if nargin < 7 || isempty(dryrun)
    dryrun = false;
end


 gettank = invoke(TT,'Opentank',fullfile(tpath,tankname),'r');
if gettank == 0
    pth = fileparts(fullfile(tpath,tankname));
    pth2 = fileparts(pth(1:end-1));
    gettank = invoke(TT,'Opentank',pth2,'r');
    if gettank == 0
        error('Failed to open tank.')
    end
end

if nargin < 6
    FS=[];
end

if isstr(blocknames)
    blocknames = {blocknames};
end

i = 0;
if nargin < 2 || isempty(blocknames)
    while 1
       i = i+1; 
        blks{i} = invoke(TT,'QueryBlockname',i);
        if isempty(blks{i})
            blks(i) = [];
            break
        end
     
    end
    
    if length(blks) > 1 && interactive
        
    
           if use_listdlg && usejava('jvm')
                inp =listdlg('PromptString','Which block?',...
                        'ListString',blks);
            else
                inpstr = sprintf('%i. %%s\n',1:length(blks));
                inp = input(sprintf(['Which Block?\n',inpstr,':'],blks{:}));
           end
        blocknames = blks(inp);
    else
        blocknames = blks(1);
    end
else
     if isfloat(blocknames)
        for q = 1:length(blocknames)
           blockn{q} = invoke(TT,'QueryBlockName',blocknames(q));
           if isempty(blockn)
               error('There aren''t that many blocks.')
           end
        end
        blocknames = blockn;
    end
    
    
end

if nargin < 3 
    DataCode = '';
end

if nargin < 4
    chnum = [];
end

%increase maximum memory
TT.SetGlobalV('WavesMemLimit', 1024^3)

for blk = 1:length(blocknames)

    blockname = blocknames{blk};
    invoke(TT,'selectblock',blockname);

    if  isempty(DataCode)

        evntcodes = invoke(TT,'GetEventCodes',0);
        if ~isnan(evntcodes)
        for i = 1:length(evntcodes);

            evnt{i} = invoke(TT,'CodeToString',evntcodes(i));
            %if ismember(evnt{1})

        end
        else
            evnt = {'xxxx'};
        end

     %    evnt(ismember(evnt,'TICK')) = [];

        if length(evnt) > 1 && interactive

             if use_listdlg && usejava('jvm')
                    inp =listdlg('PromptString','Which Record?',...
                                'ListString',evnt);
             else
                inpstr = sprintf('%i. %%s\n',1:length(evnt));
                inp = input(sprintf(['\nWhich Record?\n',inpstr,':'],evnt{:}));

            end
            DataCode = evnt(inp);
        else
            DataCode = evnt(1);
        end
    end

    if isempty(chnum) && interactive
       if use_listdlg && usejava('jvm')
           inp = inputdlg('Which Channels? ');
            chnum = str2num(inp{1});            %#ok<ST2NM>
       else           
            chnum = input('\nWhich Channels? ');
       end
    end


    if ~iscell(DataCode)
        DataCode = {DataCode};
    end


    invoke(TT,'CreateEpocIndexing');
    for i= 0:50 
        epocCodes{i+1} = invoke(TT,'GetEpocCode',i);
        if isempty(epocCodes{i+1})
            epocCodes(i+1) = [];
            break
         end

    end
    isEventCode = ismember(DataCode,epocCodes);

    %invoke(TT,'ResetFilters');
    %invoke(TT,'SetGlobalStringV','Options','FILTERED');
    %invoke(TT,'CreateEpocIndexing');
    %invoke(TT,'setGlobals','Options=FILTERED');



    for dc = 1:length(DataCode)
     if ~isEventCode( dc)   
          if chnum == 0

            invoke(TT,'SetGlobalV','Channel',chnum);
            %q = invoke(TT,'ReadWavesV',WaveName);
            N=invoke(TT,'ReadEventsV',1e6,DataCode{dc},chnum,0,dryrun,2*dryrun,'ALL');
            if N==1e10, warning('Read Maximum Events'); end
            if N<0, warning('Failed to extract data. Most likely ran of memory error.'); end

            ch = invoke(TT, 'ParseEvInfoV', 0, N, 4);
    %         q=invoke(TT,'ParseEvV',0,N);
    %         time = invoke(TT, 'ParseEvInfoV', 0, N, 6);
    %         fs = invoke(TT, 'ParseEvInfoV', 0, N, 9);
    %         varargout{blk}(dc).dat = q;
    %         varargout{blk}(dc).t = time;
    %         varargout{blk}(dc).fs = fs;
    %         varargout{blk}(dc).chan = ch;
    %         varargout{blk}(dc).tank = tankname;
    %         varargout{blk}(dc).block = blockname;
    %         varargout{blk}(dc).code = DataCode{dc};
            chnum = unique(ch);
          end

        for i = 1:length(chnum)
            invoke(TT,'SetGlobalV','Channel',chnum(i));

            %    q(:,i) = invoke(TT,'ReadWavesV',WaveName);
    %         N=invoke(TT,'ReadEventsV',1e6,DataCode{dc},chnum(i),0,0,0,'ALL')
            chunkN = 2e4;
            N = Inf;
            q = [];
            fs = [];
            gett=0;
            while N>=chunkN
                if ~isempty(FS)
                    invoke(TT,'SetGlobalV','WaveSF',FS);
                   q = invoke(TT,'ReadWavesV',DataCode{dc});
                   N= 1;
                else
                    N=invoke(TT,'ReadEventsV',chunkN,DataCode{dc},chnum(i),0,dryrun +gett ,2*dryrun,'ALL')
                   qq =invoke(TT,'ParseEvV',0,N);
                end
                ch = invoke(TT, 'ParseEvInfoV', 0, N, 4);
                time = invoke(TT, 'ParseEvInfoV', 0, N, 6);
                fs = [fs,invoke(TT, 'ParseEvInfoV', 0, N, 9)];
                q = [q;qq(:)];
                gett = (length(q)+1)./fs(end);
            end
            if max(fs)>0
                varargout{blk}(dc).dat(:,i) = q;
                varargout{blk}(dc).t = time; %#ok<*AGROW>
                varargout{blk}(dc).fs = fs;
                varargout{blk}(dc).chan(i) = ch(1);
                varargout{blk}(dc).tank = tankname;
                varargout{blk}(dc).block = blockname;
                varargout{blk}(dc).code = DataCode{dc};
            else %%% Assume that if sampling frequency is 0, this was an event code missed by GetEpocCode.
                 warning('Detected fs=0, assuming store is an event store!');
                 %%% Re-read from channel 0
                 N=invoke(TT,'ReadEventsV',1e6,DataCode{dc},0,0,dryrun,2*dryrun,'ALL')
                 q =invoke(TT,'ParseEvV',0,N);
                  time = invoke(TT, 'ParseEvInfoV', 0, N, 6);
                  varargout{blk}(dc).evnt = q;
                varargout{blk}(dc).time = time;
                varargout{blk}(dc).tank = tankname;
                varargout{blk}(dc).block = blockname;
                varargout{blk}(dc).code = DataCode{dc};                
            end
        end
         %end
     else
        chunkN = 5e4;
        N = Inf;
        w = [];
        gett = 0;
        while N>=chunkN
            ww = invoke(TT,'GetEpocsV',DataCode{dc},gett,0,chunkN);   
            gett = ww(2,end)+1e-9;
        %    w = invoke(TT,'GetEpocsV','Evnt',0,0,1e6);   
        %    w = invoke(TT,'GetEpocsV','EVNT',0,0,1e6);
            w = [ww,w];
            N = size(w,2);
        end
        varargout{blk}(dc).evnt = w(1,:);
        varargout{blk}(dc).time = w(2,:);
        varargout{blk}(dc).tank = tankname;
        varargout{blk}(dc).block = blockname;
        varargout{blk}(dc).code = DataCode{dc};
     end    

    end

end

if nargout < length(blocknames)+1
 invoke(TT,'CloseTank');
 invoke(TT,'ReleaseServer');
else
    varargout{end+1} = TT;
end
if nargout > length(blocknames)+1
    fid = fopen([mfilename,'.m']);
    varargout{end+1} = fread(fid,'uchar=>char')';
    fclose(fid);
end
% 
% 
%  

function T = tdtTanks(TT)
%
%Generates a list of registered tdt tanks

%TT=actxcontrol('TTank.x');
%invoke(TT,'ConnectServer','Local','Me');

i = 0;
T = {''};
while 1
    
    tname = invoke(TT,'GetEnumTank',i);
    if isempty(tname)
        break
    end
    T{i+1} = tname;
    i = i+1;
end


%Generates a list of data tanks