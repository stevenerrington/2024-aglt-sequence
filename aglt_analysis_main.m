%% Configuration
% Setup workspace

%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
close all; clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
[ephysLog, stimulusLog, spike_log] = import_exp_map();
ephysLog = clean_exp_map(ephysLog);
spike_log = clean_spike_map(spike_log);
load('session_audio_latency.mat')

% Set parameters
ops.timewin = -1000:5000;
ops.bl_win = -150:-50;
ops.freq = [2 200];
ops.sdf_filter = 'Gauss';
ops.sig_threshold = 2;
ops.min_sig_time = 50;
ops.freq = [2 60];
ops.glm_timewin = -150:10:750; % Define the time window for plotting beta weights

% Set parameters for sound specific stimuli
ops.sound_sdf_window = -200:800;
sound_onset_ms = [0, 563, 1126, 1689, 2252];
zero_offset = abs(ops.timewin(1));

% Set colors
color_pal.auditory_clu = [46 17 45; 84 0 50; 130 3 51; 201 40 62]./255;
color_pal.frontal_clu = [15 45 64; 25 71 89; 41 107 115; 62 140 132]./255;

%% Recording information
% Define indices for neurons in specific brain regions
auditory_neuron_idx = find(ismember(spike_log.area, 'auditory'));
frontal_neuron_idx = find(ismember(spike_log.area, 'frontal'));

disp([sum(strcmp(spike_log.monkey(auditory_neuron_idx),'troy')), sum(strcmp(spike_log.monkey(auditory_neuron_idx),'walt'))]);
disp([sum(strcmp(spike_log.monkey(frontal_neuron_idx),'troy')), sum(strcmp(spike_log.monkey(frontal_neuron_idx),'walt'))]);

%% Extract element aligned data
% This code snippet is designed to manage the extraction and loading of element-aligned 
% spiking data data. The script first checks if a specific file, |sdf_soundAlign_data.mat|, 
% already exists in a designated directory. If the file does not exist, the script 
% triggers the extraction process by calling a function or script named |get_sound_aligned_data|, 
% which generates and saves the SDF data. If the file is found, it simply loads 
% the pre-existing data.

if ~exist(fullfile(dirs.root,'data','sound_align','sdf_soundAlign_data.mat'))
    fprintf('Extracting sound aligned data \n')
    get_sound_aligned_data
    save(fullfile(dirs.root,'data','sound_align'),'normFR_in','-v7.3')
else
    fprintf('Loading sound aligned data \n')
    load(fullfile(dirs.root,'data','sound_align','sdf_soundAlign_data.mat'));
end


%% Population analysis
% This code snippet is designed to perform population-level analyses on neural 
% data. The first script, |pca_seq_main|, executes a principal component analysis 
% (PCA) on sequential neural recordings to reduce dimensionality and identify 
% dominant patterns in the population activity. 

pca_seq_main

% The second script, |pca_lda_id_position|, combines PCA with linear discriminant 
% analysis (LDA) to classify neural responses according to position-related 
% information. It leverages the reduced-dimensionality data from the PCA step 
% to improve discriminative power while minimizing noise.

pca_lda_id_position

%% Single-unit analysis
% This code snippet is designed to perform detailed analyses on single-neuron 
% activity. The analyses are divided into three major components:
%
% 1. GLM analysis and clustering:
%    - |glm_singleunit_analysis| applies a Generalized Linear Model (GLM) to 
%      each neuron's spiking activity to quantify the relationship between 
%      neural firing and task-related predictors.

glm_singleunit_analysis

%- |glm_clustering| groups neurons based on their GLM-derived response 
%      patterns, enabling identification of functionally similar cell types.

glm_clustering

% 2. Sequence and periodicity analysis:
%    - |seq_autocorr| computes the autocorrelation of spike trains to assess 
%      sequential firing patterns and rhythmicity at the single-neuron level.

seq_autocorr

% 3. Intrinsic timescales:
%    - |acf_intrinsic_timescales| estimates the decay of autocorrelations over 
%      time to quantify the intrinsic temporal dynamics of individual neurons.

acf_intrinsic_timescales
