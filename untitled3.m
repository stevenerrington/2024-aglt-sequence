clear all; clc

load('/Users/stevenerrington/Desktop/troy-agl_t-2021-10-08.mat')



% Assume your data is stored in a variable called 'lfp_data' of size (channels x time x trial)
lfp_data = lfp;

% Create the FieldTrip data structure
data = [];
data.label = arrayfun(@(x) ['chan' num2str(x)], 1:size(lfp_data, 1), 'UniformOutput', false)'; % Create channel labels
data.fsample = 1000; % Sampling frequency in Hz
data.time = arrayfun(@(x) (1:size(lfp_data, 2)) / data.fsample, 1:size(lfp_data, 3), 'UniformOutput', false); % Time axis
data.trial = arrayfun(@(x) squeeze(lfp_data(:, :, x)), 1:size(lfp_data, 3), 'UniformOutput', false); % Trials

% Add additional information if needed (e.g., trialinfo)
data.trialinfo = ones(size(lfp_data, 3), 1); % Example trial info, all ones

% Check the resulting data structure
disp(data);









cfg = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.foi        = 1:2:30;
cfg.t_ftimwin  = 5./cfg.foi;
cfg.tapsmofrq  = 0.4 *cfg.foi;
cfg.toi        = -0.5:0.05:1.5;
TFRmult = ft_freqanalysis(cfg, data);

Copycfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.zlim         = [-2e-27 2e-27];
cfg.showlabels   = 'yes';
cfg.layout       = 'CTF151_helmet.mat';
cfg.colorbar     = 'yes';
figure
ft_multiplotTFR(cfg, TFRmult)
