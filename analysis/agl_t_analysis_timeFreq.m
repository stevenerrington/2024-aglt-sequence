
%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
[ephysLog, stimulusLog, spike_log] = import_exp_map();
ephysLog = clean_exp_map(ephysLog);

%%
session_i = 40;
load(fullfile(dirs.mat_data,ephysLog.session{session_i}))

%% Analysis
aligntime = event_table.stimulusOnset_ms;

% Setup LFP data for the time frequency analysis
ops.timewin = [-1000:5000];
ops.freq = [1 100];
ops.ch_extract = [1:32];
[lfp_aligned, lfp_array] = get_lfp_aligned(lfp,aligntime,ops);

% Run time-frequency power extraction limited to trials of interest
ops.tf_trials = find(~isnan(event_table.rewardOnset_ms) & strcmp(event_table.cond_label,'viol'));
tf_data = get_timefrequency(lfp_aligned, ops);

in_lfp = squeeze(lfp_array(1,:,:))';
fLFP = lfp_filter(nanmean(in_lfp), 2 , 8 ,1000)


figuren; hold on
plot(ops.timewin, fLFP)
%% Figure

figuren('Renderer', 'painters', 'Position', [100 100 900 400]);
ax1 = nsubplot(1,1,1,1);
contourf(ops.timewin,tf_data.frequencies,nanmean(tf_data.power,3),40,'linecolor','none');
% set(gca,'ytick',round(logspace(log10(tf_data.frequencies(1)),log10(tf_data.frequencies(end)),10)*100)/100,'yscale','log','xlim',[min(ops.timewin) max(ops.timewin)])

colorscale = flipud(cbrewer('div','RdBu',100));
colorscale(colorscale<0) = 0;
colormap(colorscale)
colorbar
