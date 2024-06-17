% use this to epoch Christopher's HOSD data in the same way as done for
% MountainSort (in order to make plotting easier)
clear all
%% inputs
binSize = 0.05; % specify bin size for PSTH (in seconds)
plotInHz = 1; % boolean for plotting PSTH in Hz
combinedFlag = 1; % boolean to plot all combined multi-units for each channel
maxEventsPerPlot = 6;
baselineTime = [-0.5 -0.1]; % use this window (in -/+ s) to determine auditory responsiveness - will be >= 1st & < 2nd
plotSDbothways = 0; % boolean to plot lower standard deviation threshold
audRespThresh = 3; % threshold (in standard deviations) of auditory responsiveness
rmNoisy = 0; % boolean to remove noisy trials
noisyTThresh = 2; % threshold for noisy trial removal (in SD)
markTimes = [0]; % mark times to plot
% use this to plot a select time window
subEpoch = [];
chanName = {'chan001','chan002','chan003'};
blockName = {'274-002'};
timeStampMat = {'274_002'}; % corresponding to above cell array
evntNames = {'Voice On','Voice Off','Limb On','Limb Off'}; % corresponding to exact number and order of events

for eachBlock = 1:numel(blockName)
    for eachChan = 1:numel(chanName)
        % christopher data folder
        HOSDfolder = ['\\172.30.13.77\HBRL_Upload\for_JBerger\JBData\STN\SpikeData\hosd_sorting\' blockName{eachBlock} '_ch001-003\'];
        % load timestamps
        load(['\\172.30.13.77\HBRL_Upload\for_JBerger\JBData\STN\SpikeData\Timestamps\' timeStampMat{eachBlock} '.mat'],'FIDX')
        % load data
        load([HOSDfolder 'HOSD_spike_sort_' chanName{eachChan} '.mat'])
        
        %% pull in triggers for this block and sort
        preStimSecs = 0.5;
        postStimSecs = FIDX.postStimSecs;
        events = unique(FIDX.evnt);
        eventTimes = cell(1,numel(events));
        evnt1 = find(FIDX.evnt==1);
        evnt3 = find(FIDX.evnt==3);
        evnt1EndTimes = FIDX.trialDur(1:length(evnt1));
        evnt3EndTimes = FIDX.trialDur(length(evnt1)+1:length(evnt1)+length(evnt3));
        for ii = 1:numel(events)
            eventTimes{ii} = FIDX.time(FIDX.evnt==events(ii));
        end
        
        % find unique waveform cluster numbers
        waveClusNums = unique(waveform_cluster);
        numEvents = numel(events);
        for i = 1:numel(waveClusNums)
            allSpikes = [];
            figsReq = ceil(numel(events)/maxEventsPerPlot);
            for figSelect = 1:figsReq
                waveF(figSelect) = figure;
                set(waveF(figSelect),'units','normalized','outerposition',[0 0 1 1])
            end
            
            clf
            figNum = 1;
            subplotNum = 1;
            
            indxClus = find(waveform_cluster==waveClusNums(i));
            spikeTimesClus = spike_times(indxClus);
            for ii = 1:numEvents
                timeStart = eventTimes{ii} - preStimSecs;
                timeEnd = eventTimes{ii} + postStimSecs; % in seconds, after event
                eachevntSpike = [];
                singleevntSpike = [];
                numTperEvent(ii) = length(timeStart);
                for epochi = 1:length((timeStart))
                    sp = timeStart(epochi);
                    ep = timeEnd(epochi);
                    tmp = spikeTimesClus(:,spikeTimesClus>=sp & spikeTimesClus <=ep);
                    
                    % make spike times peri-event
                    tmp = tmp-(timeStart(epochi)+preStimSecs);
                    eachevntSpike{epochi} = tmp;
                end
                % store all spikes for each unique event
                allSpikes{ii} = eachevntSpike;
                
            end
            currBlockSpikes = allSpikes;
            
            %% Calculate bins
            durEpoch = postStimSecs+preStimSecs;
            numHistC = durEpoch/binSize;
            histEdges = linspace(-preStimSecs,postStimSecs,numHistC);
            
            %% plotting
            
            for ii = 1:numEvents
                
                spikeCounts = [];
                edges = [];
                spikeTimes = {};
                %% UPTO HERE
                % calculate PSTH per trial (and then average)
                for numTrials = 1:numel(currBlockSpikes{ii})
                    allSpikes = [];
                    allSpikes = currBlockSpikes{ii}{numTrials}; % for current event
                    spikeTimes{numTrials} = allSpikes;
                    [spikeCounts(numTrials,:), edges(numTrials,:)] = histcounts(allSpikes,histEdges);
                end
                set(0, 'CurrentFigure', waveF(figNum))
                if subplotNum > maxEventsPerPlot
                    % set all axes the same
                    linkaxes(ax(:),'xy')
                    linkaxes(rastx(:),'xy')
                    linkaxes([ax(:) rastx(:)],'x')
                    figNum = figNum+1;
                    subplotNum = 1;
                end
                rastx(subplotNum,:) = subplot(2,numEvents,subplotNum);
                
                % plot end times if event == 1 or 3
                if ii == 1
                    % sort according to event Times
                    [evnt1EndTimes_sorted, evnt1EndTimes_order] = sort(evnt1EndTimes);
                    spikeTimes = spikeTimes(evnt1EndTimes_order);
                    plotSpikeRaster(spikeTimes,'PlotType','vertline');
                    hold on
                    scatter(evnt1EndTimes_sorted,1:length(evnt1EndTimes_sorted),'r','filled');
                    ylabel('Trial #')
                elseif ii == 3
                    % sort according to event Times
                    [evnt3EndTimes_sorted, evnt3EndTimes_order] = sort(evnt3EndTimes);
                    spikeTimes = spikeTimes(evnt3EndTimes_order);
                    plotSpikeRaster(spikeTimes,'PlotType','vertline');
                    hold on
                    scatter(evnt3EndTimes_sorted,1:length(evnt3EndTimes_sorted),'r','filled');
                else
                    plotSpikeRaster(spikeTimes,'PlotType','vertline');
                    hold on
                end
                
                if exist('evntNames','var')
                    title(['Event ID:' evntNames(ii)])
                else
                    title(['Event ID:' num2str(events(ii))])
                end
                hsp1 = get(gca, 'Position');                     % Get 'Position' for raster
                
                %% Remove particularly noisy trials if desired
                if rmNoisy == 1
                    noiseThresh = mean(spikeCounts,'all')+(noisyTThresh*std(spikeCounts,[],'all'));
                    rmTrials = mean(spikeCounts,2)>noiseThresh;
                    spikeCounts = spikeCounts(~rmTrials,:);
                    disp([num2str(sum(rmTrials)) ' trial(s) removed'])
                end
                
                ax(subplotNum,:) = subplot(2,numEvents,subplotNum+numEvents);
                
                summedspikes = sum(spikeCounts,1);
                edges = mean(edges,1);
                %% Convert to Hz if desired
                if plotInHz == 1
                    spikeCounts = (summedspikes./binSize)/size(spikeCounts,1);
                    histogram('BinEdges', edges, 'BinCounts',spikeCounts,'FaceColor','k')
                    ylabel('Mean firing rate (Hz)')
                else
                    spikeCounts = bar(summedspikes./size(spikeCounts,1));
                    ylabel(['Average spikes per ' num2str(binSize) ' secs'])
                end
                
                %% Make plot prettier
                xlabel('Time (s)')
                ylims = get(gca,'YLim');
                hold on
                if ~isempty(subEpoch)
                    xlim(subEpoch)
                end
                hsp2 = get(gca, 'Position');                     % Get 'Position' for (2,1,2)
                set(gca, 'Position', [hsp2(1:2)  hsp1(3:4)]);    % Use 2*(2,1,1) Height for (2,1,2)
                
                % plot triggers
                for trigPlot = 1:numel(markTimes)
                    plot([markTimes(trigPlot) markTimes(trigPlot)],ylims,'--r')
                end
                %% Determine auditory responsiveness and plot threshold
                baseFiring = spikeCounts(edges>=baselineTime(1) & edges<baselineTime(2));
                responsiveThresh = mean(baseFiring)+(audRespThresh*std(baseFiring));
                xlims = get(gca,'XLim');
                plot(xlims,[responsiveThresh responsiveThresh],'--g')
                if plotSDbothways == 1
                    plot(xlims,[mean(baseFiring)-(audRespThresh*std(baseFiring)),...
                        mean(baseFiring)-(audRespThresh*std(baseFiring))],'--c')
                end
                % increment subplot number
                subplotNum = subplotNum+1;
            end
            
            % set all axes the same
            linkaxes(ax(:),'xy')
            linkaxes(rastx(:),'xy')
            linkaxes([ax(:) rastx(:)],'x')
            
            for figSelect = 1:figsReq
                img = getframe(waveF(figSelect));
                imwrite(img.cdata, [HOSDfolder,...
                    blockName{eachBlock} '_' chanName{eachChan} '_clus_' num2str(i)  '_fig_' num2str(figSelect) '.png'], 'png');
            end
            
        end
    end
end