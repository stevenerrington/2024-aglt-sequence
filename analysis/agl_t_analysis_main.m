%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- AGLt analysis main  -------------------------------------------------
      S P Errington, 2024

This script houses analysis pertaining to the analysis of pre-processed
AGLt data.

///////////////////////////////////////////////////////////////////////////
%} 

%% Setup workspace
% Load in matlab data
datafile = 'troy-agl_t-2021-07-21.mat';
load(fullfile(dirs.mat_data,datafile))

% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patch faulty channel
record_idx = find(strcmp(ephysLog.session,datafile(1:end-4)),1);
fault_ch_idx = ephysLog.faulty_ch(record_idx);
lfp = patch_fault_ch(lfp,fault_ch_idx);

%%  Local field potential analyses
% Time frequency & ERSP -----------------------
agl_t_analysis_ERSP








%% Laminar analyses

laminar_info.auditory = get_laminar_info([1:16], lfp_aligned, event_table);
laminar_info.vlpfc = get_laminar_info([17:32], lfp_aligned, event_table);

get_laminar_plotSummary