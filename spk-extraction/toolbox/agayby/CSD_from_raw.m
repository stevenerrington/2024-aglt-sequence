%%20210409 Finished first version of the script
clc
clear all
%to make file selection with a gui
[file,path,indx] = uigetfile('*.ncs','Select NCS file');

%For more script function
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


%%
tic
% parts = strsplit(path, '\');
% DirPart = parts{end-1};
[~,fileNew] = fileparts(file); 
% fileNew = fullfile(path,[fileNew '.mat']); %% end-4 to remove the extension
% save(fileNew,'data','event','hdr')
fprintf('\nOutput saved, took %2.2f seconds\n', toc);




for ind = 1:length(event)
    event(ind).timestampSeconds = (event(ind).timestamp - double(hdr.FirstTimeStamp)) *1e-6;
end

time = data.time{1,1}(end)
time_event = event(end).timestampSeconds
%% using marker from the serial input
temp = find([event.port]==2);

temp = event(temp);

marker = dec2bin([temp.ttls]);
marker = bin2dec(marker(:,2:8));
onset = find(marker==5);
trl = [temp(onset).timestampSeconds]';
trl = round(trl * hdr.TimeStampPerSample);

trl = [trl-1000 trl -500*ones(size(trl))];
%% using sound input
% temp = find([event.port]==0);
% sound = event(temp);
% temp = diff([sound.timestampSeconds]);
% temp = find(temp>0.8);
% trl = [sound(temp).timestampSeconds];
% 
% trl = round(trl * hdr.TimeStampPerSample)';
% 
% trl = [trl-1000 trl -500*ones(size(trl))];
%%
cfg = [];
cfg.trl = trl;
lfp = ft_redefinetrial(cfg,data);

%%
% cfg=  [];
% figure(1)
% cfg.layout = 'vertical';
% cfg.xlim = [-0.2 inf];
% % cfg.trials = [1:5];
% 
% subplot(2,2,1)
% % ft_multiplotER(cfg,lfp)
% ft_singleplotER(cfg,lfp)
% xlabel('Time (s)')
% ylabel('Voltage (uV)')
% grid on
% hold on
% xline(0)
% legend('LFP','Sound Onset')
%%
% 
% cfg_fft = [];
% cfg_fft.method = 'mtmfft';
% cfg_fft.foilim = [0 80];
% cfg_fft.taper = 'dpss';
% cfg_fft.tapsmofrq = 4;
% cfg_fft.trials = [1:5];
% 
% freq = ft_freqanalysis(cfg_fft,lfp);
% 
% cfg = [];
% cfg.layout = 'vertical';
% figure(2)
% subplot(2,1,1)
% ft_singleplotER(cfg,freq)
% grid on
% % title('R1 Ref 02 as ref')
% 
% %%
%     cfg_tf = [];
%     cfg_tf.method = 'wavelet';
%     cfg_tf.width = 7;
%     cfg_tf.foi = [1:1:50];
% %     cfg_tf.toi = [0:0.01:30];
%     cfg_tf.toi = 'all';
%     tfreq = ft_freqanalysis(cfg_tf,lfp);
%     %%
%     cfg = [];
% %     cfg.layout = 'vertical';
%     cfg.ylim = [0 50];
%     figure(1)
%     subplot(2,2,3)
%        colormap(jet)
%     ft_singleplotTFR(cfg,tfreq)
%% CSD calculation from LFP
    cfg = []; 
%     cfg.trials = [5:20];
    lfp_csd = ft_selectdata(cfg,lfp);
    
    [csd,icsd,time] = calculateCSD_nylx(lfp_csd);
%%
    figure(1)
    num_ch = 16;
    ch = 0:(num_ch+1)/200:(num_ch+1)-((num_ch+1)/200);
    subplot(2,2,[2 4])
    imagesc(time,ch,1.*(icsd))

    colormap(flipud(jet))
    colorbar
    xlabel('Time (seconds)')
    yt = get(gca,'YTick');
% figure
%     saveas(gcf,[file '_csd_' my_datetimestr '.fig']) 
%     figure
%     subplot(1,2,2)
%     imagesc(time,[num_ch:1],(csd))
%        colormap(flipud(jet))
%     colorbar
%     xlabel('Time (seconds)')
%     yt = get(gca,'YTick');



% %% Plotting MUA
%     cfg = [];
%     cfg.layout = 'vertical';
%     cfg.baseline = [-0.1 0];
%     figure
%     ft_multiplotER(cfg,muax)