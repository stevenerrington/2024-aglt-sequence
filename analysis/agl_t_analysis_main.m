clear all; clc
dirs = set_directories();
datafile = 'troy-agl_t-2021-11-05.mat';
load(fullfile(dirs.mat_data,datafile))

% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find event information to align data on
aligntime = event_table.stimulusOnset_ms;

% Run alignment algorithms
ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';
[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);
lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

%%  Local field potential analyses
% Time frequency -----------------------
% This script will cycle through the LFP 
agl_t_analysis_timeFreq

% Intertrial phase clustering ----------
agl_t_analysis_intertrialphasecoherence

%% Laminar analyses

laminar_info.auditory = get_laminar_info([1:16], lfp_aligned, event_table);
laminar_info.vlpfc = get_laminar_info([17:22,24:32], lfp_aligned, event_table);

get_laminar_plotSummary