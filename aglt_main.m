
%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- Kikuchi lab main script -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%} 

%% Workspace configuration and setup
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
ephysLog = import_exp_map();
ephysLog = clean_exp_map(ephysLog);

%% Data extraction
% This script will extract data from noted recording system, namely:
% lfp data, spike times, spike waveforms, spike metrics, event times and
% trials. They will then be saved to the matlab directory defined in
% set_directories().

aglt_multisession_loop






%% Analysis
agl_t_analysis_main


%% Appendix
kikuchi_neural_HOSD % A temporary script using HOSD as an alternative to KS
laminar_summary_loop % Loop through neural data and get a laminar summary figure
kikuchi_nlx_extract      % A standalone script that was developed to extrac Neuralynx data - now in function form
kikuchi_tdt_extract      % A standalone script that was developed to extrac TDT data - now in function form
