
%% 
%{ 
----- Kikuchi lab main script -----------------------------------------
      S P Errington, 2024
%} 

%% Workspace configuration and setup
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();


%% Data extraction
% This script will extract data from noted recording system, namely:
% lfp data, spike times, spike waveforms, spike metrics, event times and
% trials. They will then be saved to the matlab directory defined in
% set_directories().

kikuchi_nlx_extract % Neuralynx
% ---------         % TDT
% ---------         % NeuroNexus


kikuchi_neural_HOSD % A temporary script using HOSD as an alternative to KS


%% Analysis
agl_t_analysis_main

