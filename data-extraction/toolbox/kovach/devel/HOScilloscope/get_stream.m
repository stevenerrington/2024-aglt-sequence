

serverName = '172.28.253.13'

succeeded = NlxConnectToServer(serverName);
if succeeded ~= 1
    disp(sprintf('FAILED connect to %s. Exiting script.', serverName));
    return;
else
    disp(sprintf('Connect successful.'));
end
%[succeeded, reply] = NlxSendCommand('-StartAcquisition')
serverIP = NlxGetServerIPAddress();
disp(sprintf('Connected to IP address: %s', serverIP));

serverPCName = NlxGetServerPCName();
disp(sprintf('Connected to PC named: %s', serverPCName));

serverApplicationName = NlxGetServerApplicationName();
disp(sprintf('Connected to the NetCom server application: %s', serverApplicationName));


%get a list of all objects in Cheetah, along with their types.
[succeeded, cheetahObjects, cheetahTypes] = NlxGetDASObjectsAndTypes;
if succeeded == 0
    disp 'FAILED get cheetah objects and types'
else
    disp 'PASSED get cheetah objects and types'
end
%%

getobjects = listdlg('PromptString','Select Objects','ListString',cheetahObjects,'SelectionMode','multi');

%open up a stream for all objects
for index = getobjects
    succeeded = NlxOpenStream(cheetahObjects(index));
    if succeeded == 0
        disp(sprintf('FAILED to open stream for %s', char(cheetahObjects(index))));
        break;
    end
end;
if succeeded == 1
    disp 'PASSED open stream for all current objects'
end

%send out an event so that there is something in the event buffer when
%this script queries the event buffer  You can use NlxSendCommand to send
%any Cheetah command to Cheetah.
[succeeded, cheetahReply] = NlxSendCommand('-PostEvent "HOScilloscope connected" 10 11');
if succeeded == 0
    disp 'FAILED to send command'
else
    disp 'PASSED send command'
end

nbin = 30;

%%
heatbin = zeros(length(getobjects),nbin);
fig = figure;
lowpass = 150; 
windur = 1;
stepsize = .1;
createnew=true;
scrollbufferN=2;
ncascade =3;
if createnew
    clear hos
end
for k = 1:length(getobjects)
    getObject = cheetahObjects{getobjects(k)};
    [succeeded,dataArray, timeStampArray, channelNumberArray, samplingFreqArray, numValidSamplesArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewCSCData(getObject);
    t0 = timeStampArray(1);
    newdata = double(dataArray(:));
    scale = std(newdata);

    sampling_rate = double(samplingFreqArray(1));
    N = round(sampling_rate*windur);
   
    if succeeded && createnew
         for kk = 1:ncascade
   
            hos(k,kk) = hosobject(3,N,sampling_rate,lowpass);
            hos(k,kk).get_input(newdata);
            hos(k,kk).hos_learning_rate = .0001;
            hos(k,kk).filter_adaptation_rate = .001;
            hos(k,kk).burnin = 1/hos(k).hos_learning_rate;
         end
            [xr,xf] = hos(k,1).reconstruct(newdata);
    elseif ~succeeded
        error('Failed for %s',cheetahObjects{getobjects(k)})
    end
    wf = fftshift(hos(k).freqindx.Bfreqs{1});
    handles.ax(k,1) = subplot(length(getobjects),3,3*(k-1)+1);
    handles.im(k,1) = imagesc(wf,wf,fftshift(abs(hos(k).bicoh)));
    tth(k) = title(sprintf('EDF: %0.2f',hos(k).EDF));
    axis image
    colorbar
    handles.ax(k,2) = subplot(length(getobjects),3,3*(k-1)+2);
    title(cheetahObjects{ getobjects(k)});
    handles.pl1(k,1:ncascade+1) = plot(ifftshift(hos(k).sampt)/sampling_rate,[fftshift([hos(k,:).waveform],1),fftshift(hos(k).shiftbuffer)]);
    set(handles.pl1(k,1),'color',[1 1 1]*.75)
    axis tight
    ylim([-1 1]*3*scale)
    title(cheetahObjects{ getobjects(k)});
     handles.ax(k,3) = subplot(length(getobjects),3,3*(k-1)+3);
    
     scrollbuffer{k} = zeros(N*scrollbufferN,1);
      t = (0:length( scrollbuffer{k})-1)./hos(k).sampling_rate+double(timeStampArray(1)-t0)/1e6;
    
    [xr,xf] = hos(k).reconstruct(scrollbuffer{k});
%      handles.pl2(k,1+(1:1+ncascade)) = plot(t,[scrollbuffer{k},xr(:,ones(1,ncascade))]);
    handles.pl2(k,1:2) = plot(t,[scrollbuffer{k},xr]);
    set(handles.pl2(1),'color',[1 1 1]*.5);
    ylim([-1 1]*6)
   heatbin(k,end) = skewness(xf);
   axis tight
end
if createnew
    for k = 1:length(getobjects)
        for kk=1:ncascade
        hos(k,kk).reset;
%          hos(k,kk).highpass = 100;
        end
    end
end

% subplot(length(getobjects),3,3)
% imhb = imagesc(heatbin);
% numberOfPasses = 1e3;
% colorbar
%%
% newdatas = repmat({[]},1,length(getobjects));
nhist = 500;
xhist = zeros(length(scrollbuffer{1}),nhist,length(hos));
stepn = zeros(size(hos));
numberOfPasses = 1e6; 
for pass = 1:numberOfPasses
   tic;
    for k = 1:length(getobjects)
        objectIndex = getobjects(k);
        objectToRetrieve = char(cheetahObjects(objectIndex));
        %determine the type of acquisition entity we are currently indexed
        %to and call the appropriate function for that type
        if strcmp('CscAcqEnt', char(cheetahTypes(objectIndex))) == 1
            [succeeded,dataArray, timeStampArray, channelNumberArray, samplingFreqArray, numValidSamplesArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewCSCData(objectToRetrieve);
        end
%          newdata = zscore(double(dataArray(:)));
           newdata = double(dataArray(:));
         
%          newdatas{k} = [newdatas{k};newdata];
        ndl = length(newdata);
%                 hos(k).get_input(newdata);

        if ndl < scrollbufferN*hos(k).buffersize
            scrollbuffer{k}(1:end-ndl,1)=  scrollbuffer{k}(ndl+1:end);
            scrollbuffer{k}(end-ndl+1:end,1) = newdata;
        else
            scrollbuffer{k} = newdata(1:scrollbufferN*hos(k).buffersize);
        end
%         if succeeded 
%              hos(k).get_input(newdata);
%          else
%             warning('Failed for %s',cheetahObjects{getobjects(k)})
%         end
        stepn(k)=stepn(k)+length(newdata);
        if  stepn(k) >= hos(k).buffersize*stepsize
            indata = scrollbuffer{k};
                set(handles.im(k,1), 'cdata',fftshift(abs(hos(k,1).bicoh)));
                set(tth(k),'string',sprintf('EDF: %0.2f',hos(k,1).EDF))
                set(handles.pl1(k,1), 'ydata',fftshift(hos(k,1).shiftbuffer));
               t = (0:length(scrollbuffer{k})-1)./hos(k).sampling_rate+double(timeStampArray(1)-t0)/1e6;
               set(handles.pl2(k,1), 'xdata',t', 'ydata',indata); 
            for kk = 1:ncascade
                  hos(k,kk).get_input(indata);
                     [xr(:,kk),xf(:,kk)] = hos(k).reconstruct(indata);
                set(handles.pl1(k,2+kk-1), 'ydata',fftshift(hos(k,kk).waveform));
%                 heatbin(k,end,kk) = skewness(xf);
%                set(handles.pl2(k,2+kk), 'xdata',t','ydata',xr); 
                 indata = indata-xr(:,kk);
                xhist(:,mod(pass-1,nhist)+1,k,kk) =indata;
            end
               set(handles.pl2(k,2), 'xdata',t','ydata',sum(xr,2)); 

                xlim(handles.ax(k,3),t([1 end]))
                stepn(k)=0;
        end

    end
          drawnow
    
%     set(imhb,'cdata',heatbin)
%     heatbin(:,1:end-1) = heatbin(:,2:end);
   
   while toc<windur*stepsize
   end
%    hos(1).EDF
%    hos(1).window_number
end

%[succeeded, reply] = NlxSendCommand('-StopAcquisition')
