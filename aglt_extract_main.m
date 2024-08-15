
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
spike_log = clean_spike_map(spike_log);
load('session_audio_latency.mat')

%% Data extraction ////////////////////////////////////////////////////////
% This script will extract data from noted recording system, namely:
% lfp data, spike times, spike waveforms, spike metrics, event times and
% trials. They will then be saved to the matlab directory defined in
% set_directories().

% [! main extraction] aglt_multisession_loop

%% Analysis
agl_t_analysis_spkModulation