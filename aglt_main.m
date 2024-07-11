
%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- Kikuchi lab main script -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%} 

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

%% Data extraction ////////////////////////////////////////////////////////
% This script will extract data from noted recording system, namely:
% lfp data, spike times, spike waveforms, spike metrics, event times and
% trials. They will then be saved to the matlab directory defined in
% set_directories().

% [! main extraction] aglt_multisession_loop


%% Data curation /////////////////////////////////////////////////////////
% First, we will determine which of our sessions may have recordings that
% are across the cortical sheet.

%% Analysis //////////////////////////////////////////////////////////////
% This script will run scripts related to the analysis of AGLt data. 
% (1) Time frequency/ERSP 
agl_t_analysis_main


%% Appendix
% This section houses scripts that were used in development that may serve
% as useful references.
kikuchi_neural_HOSD % A temporary script using HOSD as an alternative to KS
laminar_summary_loop % Loop through neural data and get a laminar summary figure
kikuchi_nlx_extract      % A standalone script that was developed to extract Neuralynx data - now in function form
kikuchi_tdt_extract      % A standalone script that was developed to extract TDT data - now in function form
