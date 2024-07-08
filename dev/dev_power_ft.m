
%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
[ephysLog, stimulusLog] = import_exp_map();
ephysLog = clean_exp_map(ephysLog);

filename = 'troy-agl_t-2021-10-08';
load(fullfile('C:\KIKUCHI-LOCAL\data\ephys\mat\', [filename '.mat']))

% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patch faulty channel
record_idx = find(strcmp(ephysLog.session,filename),1);
fault_ch_idx = ephysLog.faulty_ch(record_idx);
lfp = patch_fault_ch(lfp,fault_ch_idx);

% Assume your data is stored in a variable called 'lfp_data' of size (channels x time x trial)
lfp_data_in = lfp;

% Get aligned signal
aligntime = event_table.stimulusOnset_ms;
ops.timewin = [-1000:5000];

% Theta extraction -----------------------------------------------
ops.freq = [4 8];
[theta_lfp_data] = get_lfp_aligned(lfp_data_in,aligntime,ops);

morletParameters.samplingFreq = 1000;
morletParameters.frequencies = 4:8;
morletParameters.cycle = 7;

clear morletLFP*

for ch_i = 1:32
    [morletLFP] = convMorletWaveform(theta_lfp_data.(['lfp_' int2str(ch_i)]),morletParameters);
    morletLFP = nanmean(morletLFP,3);
    morlet_out(:,:,ch_i) = morletLFP;
end

morlet_out_theta = morlet_out;


figuren; hold on

for stimuli_i = 1:16
    trials_in = []; trials_in = find(event_table.cond_value == stimuli_i);
    a = subplot(2,1,1); hold on
    morletLFP_all_1 = []; morletLFP_all_1 = nanmean(morlet_out_theta(trials_in,:,[1:16]),3);
    plot(ops.timewin, nanmean(morletLFP_all_1))
    b = subplot(2,1,2); hold on
    morletLFP_all_2 = []; morletLFP_all_2 = nanmean(morlet_out_theta(trials_in,:,[16:32]),3);
    plot(ops.timewin, nanmean(morletLFP_all_2))

end

    linkaxes([a b],'y')




% Convert data into FieldTrip format
ft_data = klm_to_ft(lfp_data, ops);

















%% Extract iSpline CSD
cfg                 = [];
cfg.vartrllength    = 1;
cfg.channel         = 1:length(data.label);
CSDdata             = ft_timelockanalysis(cfg, data);
% 
cfg = [];
cfg.baseline = [-0.50 0];
CSDdata = ft_timelockbaseline(cfg,CSDdata);
cfg = [];
cfg.latency = [-0.5 1.2];
CSDdata = ft_selectdata(cfg,CSDdata);

dt                  = 0.5;
diam                = 0.5;
cond                = 0.4;
cond_top            = 0.4;
el_d                 =  0.2;
el_pos              = (0.1:el_d:((length(data.label))*el_d))*1e-3;

gauss_sigma         = 0.1*1e-3;
filter_range        = 5*gauss_sigma;

Fcs                 = F_cubic_spline(el_pos,diam,cond,cond_top);
[zs,CSD_cs]         = make_cubic_splines(el_pos,CSDdata.avg,Fcs);
[zs,CSD_cs]         = gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
iCSD                = CSD_cs/10^6;



%% 
cfg = [];
cfg.channel    = ft_data.label;
cfg.method     = 'wavelet';
cfg.width      = 3;
cfg.output     = 'pow';
cfg.foi        = 1:1:100;
cfg.toi        = -1.0:0.05:5.0;
TFRwave = ft_freqanalysis(cfg, ft_data);

cfg              = [];
cfg.baseline     = [-0.5 -0.0];
cfg.baselinetype = 'db';
cfg.maskstyle    = 'saturation';
cfg.channel      = ft_data.label;
cfg.interactive  = 'no';
figuren
ft_singleplotTFR(cfg, TFRwave);
