clc
close all
clear all
[file,path,indx] = uigetfile('*.mat','Select mua file');

load(fullfile(path,file))


% temp= [event.timestamp] - double(hdr.FirstTimeStamp);

for ind = 1:length(event)
    event(ind).timestampSeconds = (event(ind).timestamp - double(hdr.FirstTimeStamp)) *1e-6;
end

time = data.time{1,1}(end)
time_event = event(end).timestampSeconds


%% This would be using the code to extract stimulus onset
temp = find([event.port] == 2);
trl_event = event(temp);
trl_event(1).code = [];

temp = dec2bin([trl_event.value]);
temp = bin2dec(temp(:,2:8));
% temp = num2cell(temp);
temp = find(temp == 9);
trlbegin = [trl_event(temp).timestampSeconds];
% trlbegin = trlbegin(1:2:end);
trlbegin = round( trlbegin * hdr.Fs )'; %that's the onset of stimulation. 
trlend = trlbegin + 1.5 * hdr.Fs;
trlbegin = trlbegin - 0.5 * hdr.Fs;
offset = 0.5*hdr.Fs * ones(size(trlend));

trl = [trlbegin trlend -offset];


%%
cfg = [];
cfg.trl = trl;
mua = ft_redefinetrial(cfg,data);

cfg = [];
cfg.resamplefs = 1000;
muae = ft_resampledata(cfg,mua);
% 
for ind=1:length(muae.trial)
    muae.trial{1,ind} = ft_preproc_rectify(muae.trial{1,ind});
end

%%
cfg=  [];
% figure(1)
% figure
cfg.layout = 'vertical';
% cfg.xlim = [-0.2 inf];
% cfg.trials = [1:50];
cfg.baseline = [-0.5 -0.1];
% subplot(2,2,1)
ft_multiplotER(cfg,muae)
% ft_singleplotER(cfg,lfp)
xlabel('Time (s)')
ylabel('Voltage (uV)')
grid on
% hold on
% xline(0)
% legend('LFP','Sound Onset')
figure(2)
cfg.channel = [8];
cfg.xlim = [-0.05 0.4];
ft_singleplotER(cfg,muae)
xlabel('Time (s)')
ylabel('MUAe (uV)')
%%
