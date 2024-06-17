clear all; clc
dirs = set_directories();
load(fullfile(dirs.mat_data,"troy-agl_t-2021-11-05.mat"))

%% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Align data to event
% Define parameters
ops.event_port = 2;
events = get_agl_t_trials(event_table, ops);

aligntime = events.stimulusOnset_ms;


ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);
lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

%get_laminar_profile([1:16], lfp, events)
% issues: clim just 





% Try time-frequency plot

% should be a theta? Kikuchi et al

