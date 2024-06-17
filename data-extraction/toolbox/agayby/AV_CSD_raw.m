%%20210409 Finished first version of the script
clc
clear all
close all
%to make file selection with a gui
% [file,path,indx] = uigetfile('*.ncs','Select NCS file');
[file,path,indx] = uigetfile('LFP1*.ncs','Select NCS file');
% %For more script function
% path = pwd;
% file = 'LFP1.ncs';

% tic

%% Getting individual file names for channels of interest
tic
cfg =[];
% ch = [1:16];
ch = [17:32];
temp = repmat({file(1:3)},length(ch),1);
ext = repmat({'.ncs'},length(ch),1);
ch = cellfun(@num2str, num2cell(ch), 'UniformOutput', 0)';
if isempty( strfind(file,'_') )
    cfg.fname = cellfun(@horzcat, temp,ch,ext,'UniformOutput',false);
    clearvars ch ext temp
else
    file_num = file(strfind(file,'_'):strfind(file,'.ncs')-1);
    file_num = repmat({file_num},length(ch),1);
    cfg.fname = cellfun(@horzcat, temp,ch,file_num,ext,'UniformOutput',false);
    clearvars ch ext temp file_num
end

% cfg.latency = [0 300];
%%

data = ft_read_neuralynx_interp(cfg.fname);
fprintf('\nChannel files loaded and processed, took %2.2f seconds\n', toc);
%%
tic
hdr   = ft_read_header(fullfile(path,file));

if isempty(strfind(file,'_'))    %% detecting if the file is the original one events.nev or the numbered ones events_xxxx.nev
    event_file = 'Events.nev';
else
    file_num = file(strfind(file,'_'):end-4);
    event_file = ['Events' file_num '.nev'];
end

event = ft_read_event_BA(event_file);

for i=1:length(event)
  % the first sample in the datafile is 1
  event(i).sample = (event(i).timestamp-double(hdr.FirstTimeStamp))./hdr.TimeStampPerSample + 1;
end
fprintf('\nEvents loaded and processed, took %2.2f seconds\n', toc);
%%

for ind = 1:length(event)
    event(ind).timestampSeconds = (event(ind).timestamp - double(hdr.FirstTimeStamp)) *1e-6;
end

time = data.time{1,1}(end)
time_event = event(end).timestampSeconds
%%
temp = find([event.port]==2);
% 
temp = event(temp);
marker = dec2bin([temp.ttls]);
marker = bin2dec(marker(:,2:8));

% trlbegin = find(marker ==1);
% 
% trlbegin = [temp(trlbegin).timestampSeconds]';
% trlbegin = round(trlbegin * hdr.Fs);


trlend = find(marker ==3); %3 is the reward is on
trlbegin = [];
for ind = 1:length(trlend)
    temp2 = find(marker(1:trlend(ind)) == 1);
    trlbegin(end+1) = temp2(end);
end
trialinfo = [];
for ind = 1:length(trlend)
    trltemp = marker(trlbegin(ind):trlend(ind));
    trialinfo(ind) = trltemp(find( (trltemp>50) & (trltemp < 100) ) );
end
trialinfo = trialinfo - 50;
trlbegin = [temp(trlbegin).timestampSeconds]';
trlbegin = round(trlbegin * hdr.Fs);
trlend = trlbegin + 2*hdr.Fs;

% trlend = [temp(trlend).timestampSeconds]';
% trlend = round(trlend * hdr.Fs);


% trlend = trlbegin + 500;
trl = [trlbegin trlend -0.4*hdr.Fs*ones(size(trlbegin))];
%%
cfg = [];
cfg.trl = trl;
lfp = ft_redefinetrial(cfg,data);

for ind=1:length(trl);
    lfp.trial{1,ind}(7,:) = 0.5*(lfp.trial{1,ind}(8,:)+lfp.trial{1,ind}(6,:));
end
%%
cfg=  [];
figure(1)
subplot(1,3,1)
cfg.layout = 'vertical';
% cfg.xlim = [-0.1 0.25];
% find( (trialinfo > 6) & (trialinfo <13) );
cfg.trials =find( (trialinfo <7) );;

ft_multiplotER(cfg,lfp)
subplot(1,3,2)
cfg.trials =find( (trialinfo > 6) & (trialinfo <13) );
ft_multiplotER(cfg,lfp)
subplot(1,3,3)
cfg.trials =find( trialinfo > 12 );
ft_multiplotER(cfg,lfp)
% ft_singleplotER(cfg,lfp)
xlabel('Time (s)')
ylabel('Voltage (uV)')
grid on
title('LFP')
xlabel('Time (s)')

cfg=  [];
figure(3)
subplot(1,3,1)
cfg.layout = 'vertical';
% cfg.xlim = [-0.1 0.25];
% find( (trialinfo > 6) & (trialinfo <13) );
cfg.trials =find((trialinfo > 18) & (trialinfo <=24)  );;

ft_multiplotER(cfg,lfp)
subplot(1,3,2)
cfg.trials =find( (trialinfo > 24) & (trialinfo <=30) );
ft_multiplotER(cfg,lfp)
subplot(1,3,3)
cfg.trials =find( trialinfo > 30 );
ft_multiplotER(cfg,lfp)
% ft_singleplotER(cfg,lfp)
xlabel('Time (s)')
ylabel('Voltage (uV)')
grid on
title('LFP')
xlabel('Time (s)')


%% CSD calculation from LFP
%     cfg = []; 
%     cfg.trials = find( trialinfo < 7 );
%     lfp_a = ft_selectdata(cfg,lfp);
%     
%     cfg.trials = find( (trialinfo > 6) & (trialinfo <13) );
%     lfp_v = ft_selectdata(cfg,lfp);
%     
%     cfg.trials = find( (trialinfo > 12) & (trialinfo <19) );
%     lfp_av = ft_selectdata(cfg,lfp);
%     [csd_a,icsd_a,time] = calculateCSD_nylxAV(lfp_a);
%     [csd_v,icsd_v,time] = calculateCSD_nylxAV(lfp_v);
%     [csd_av,icsd_av,time] = calculateCSD_nylxAV(lfp_av);
%% CSD calculation from LFP
    cfg = []; 
    cfg.trials = find( trialinfo < 7 );
    lfp_a = ft_selectdata(cfg,lfp);
    
    cfg.trials = find( (trialinfo > 6) & (trialinfo <13) );
    lfp_v = ft_selectdata(cfg,lfp);
    
    cfg.trials = find( (trialinfo > 12) & (trialinfo <19) );
    lfp_av = ft_selectdata(cfg,lfp);
    [csd_a,icsd_a,time] = calculateCSD_nylxAV(lfp_a);
    [csd_v,icsd_v,time] = calculateCSD_nylxAV(lfp_v);
    [csd_av,icsd_av,time] = calculateCSD_nylxAV(lfp_av);
    
    
    cfg = []; 
    cfg.trials = find((trialinfo > 18) & (trialinfo <=24) );
    lfp_a_op = ft_selectdata(cfg,lfp);
    
    cfg.trials = find( (trialinfo > 24) & (trialinfo <=30) );
    lfp_v_op = ft_selectdata(cfg,lfp);
    
    cfg.trials = find( (trialinfo > 30)  );
    lfp_av_op = ft_selectdata(cfg,lfp);
    [csd_a_op,icsd_a_op,time] = calculateCSD_nylxAV(lfp_a_op);
    [csd_v_op,icsd_v_op,time] = calculateCSD_nylxAV(lfp_v_op);
    [csd_av_op,icsd_av_op,time] = calculateCSD_nylxAV(lfp_av_op);
%%
    figure(2)
    num_ch = 16;
    ch = 0:(num_ch+1)/200:(num_ch+1)-((num_ch+1)/200);
    subplot(1,3,1)
    imagesc(time,ch,1.*(icsd_a))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD Audio')
    
    
    
    subplot(1,3,2)
    imagesc(time,ch,1.*(icsd_v))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD Visual')
    
    subplot(1,3,3)
    imagesc(time,ch,1.*(icsd_av))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD AudioVisual')
    
    
    %%
      figure(4)
    num_ch = 16;
    ch = 0:(num_ch+1)/200:(num_ch+1)-((num_ch+1)/200);
    subplot(1,3,1)
    imagesc(time,ch,1.*(icsd_a_op))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD Opto Audio')
    
    
    
    subplot(1,3,2)
    imagesc(time,ch,1.*(icsd_v_op))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD Opto Visual')
    
    subplot(1,3,3)
    imagesc(time,ch,1.*(icsd_av_op))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD Opto AudioVisual')
%%
     figure(5)
    num_ch = 16;
    ch = 0:(num_ch+1)/200:(num_ch+1)-((num_ch+1)/200);
    subplot(1,3,1)
    imagesc(time,ch,1.*(icsd_a_op - icsd_a))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD Opto Audio diff')
    
    
    
    subplot(1,3,2)
    imagesc(time,ch,1.*(icsd_v_op - icsd_v))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD Opto Visual diff')
    
    subplot(1,3,3)
    imagesc(time,ch,1.*(icsd_av_op - icsd_av))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
    title('CSD Opto AudioVisual diff')
%%
% cfg =[];
% cfg.method = 'wavelet';
% cfg.foi = [0:1:60];
% cfg.toi = [0 1.5];
% cfg.keeptrials = 'yes';
% tfreq = ft_freqanalysis(cfg,lfp);