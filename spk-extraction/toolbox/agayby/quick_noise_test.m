clear all
for ind = 1:32
    fname{ind}= ['LFP' num2str(ind) '.ncs'];
end
data = ft_read_neuralynx_interp(fname);

%%
cfg = [];
cfg.channel = [17:32];
cfg.latency = [1955 1975];
data_temp = ft_selectdata(cfg,data);
cfg_fft = [];
cfg_fft.method = 'mtmfft';
cfg_fft.foilim = [0 300];
cfg_fft.taper = 'dpss';
cfg_fft.tapsmofrq = 4;


freq = ft_freqanalysis(cfg_fft,data_temp);

    %%
    figure
    cfg = [];
    cfg.layout = 'vertical';
    cfg.xlim = [30 100];
%     cfg.channel = 1;
    ft_multiplotER(cfg,freq)