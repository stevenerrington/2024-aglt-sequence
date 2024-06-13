%%20210409 Finished first version of the script
clc
clear all
close all
%to make file selection with a gui
[file,path,indx] = uigetfile('*.ncs','Select NCS file');

% %For more script function
% path = pwd;
% file = 'LFP1.ncs';

% tic

%% Getting individual file names for channels of interest
tic
cfg =[];
ch = [1:16];
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
    trialinfo(ind) = trltemp(find( (trltemp>50) & (trltemp < 70) ) );
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
mua = ft_redefinetrial(cfg,data);
cfg = [];
cfg.resamplefs = 1000;
muae = ft_resampledata(cfg,mua);

for trl = 1:length(mua.trial)
    muae.trial{1,trl} = ft_preproc_rectify(muae.trial{1,trl});
end
%%
cfg=  [];
figure(1)
% subplot(1,2,1)
cfg.layout = 'vertical';
% cfg.xlim = [-0.1 0.25];
% cfg.trials = [1:22];
subplot(1,3,1)
cfg.trials  = find(trialinfo <7);
ft_multiplotER(cfg,muae)
title('Audio')
subplot(1,3,2)
cfg.trials  = find( (trialinfo >6) & (trialinfo <13) );
ft_multiplotER(cfg,muae)
title('Video')
subplot(1,3,3)
cfg.trials  = find(trialinfo >12);
ft_multiplotER(cfg,muae)
title('AudioVisual')
% ft_singleplotER(cfg,muae)
% xlabel('Time (s)')
% ylabel('Voltage (uV)')
% grid on
% title('LFP')
% xlabel('Time (s)')



%     yt = get(gca,'YTick');