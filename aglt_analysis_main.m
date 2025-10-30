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
% Define common parameters

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
% Count neurons

% Define indices for neurons in specific brain regions
auditory_neuron_idx = find(ismember(spike_log.area, 'auditory'));
frontal_neuron_idx = find(ismember(spike_log.area, 'frontal'));

[sum(strcmp(spike_log.monkey(auditory_neuron_idx),'troy')), sum(strcmp(spike_log.monkey(auditory_neuron_idx),'walt'))];
[sum(strcmp(spike_log.monkey(frontal_neuron_idx),'troy')), sum(strcmp(spike_log.monkey(frontal_neuron_idx),'walt'))];

% Extract & normalize sound aligned data
% This code snippet is designed to manage the extraction and loading of sound-aligned 
% spiking data data. The script first checks if a specific file, |sdf_soundAlign_data.mat|, 
% already exists in a designated directory. If the file does not exist, the script 
% triggers the extraction process by calling a function or script named |get_sound_aligned_data|, 
% which generates and saves the SDF data. If the file is found, it simply loads 
% the pre-existing data.

% get_sound_aligned_data_table

%% Extract sound aligned Spike Density Function (SDF)
if ~exist(fullfile(dirs.root,'data','sound_align','sdf_soundAlign_data.mat'))
    fprintf('Extracting sound aligned data \n')
    get_sound_aligned_data

    [normFR_in.norm_fr_soundA, normFR_in.norm_fr_soundC, normFR_in.norm_fr_soundG, normFR_in.norm_fr_soundF, normFR_in.norm_fr_soundD, normFR_in.norm_fr_soundAll] =...
    element_extract_normSoundAlign(sdf_soundAlign_data, spike_log);
    save(fullfile(dirs.root,'data','sound_align','normFR_in.mat'),'normFR_in','-v7.3')
else
    fprintf('Loading sound aligned data \n')
    load(fullfile(dirs.root,'data','sound_align','sdf_soundAlign_data.mat'));
    load(fullfile(dirs.root,'data','sound_align','normFR_in.mat'));
end

%% Population analysis
pca_seq_main
pca_lda_id_position

%% Single-unit analysis

% GLM analysis and clustering
glm_singleunit_analysis
glm_element_clustering
glm_plot_clusters

% Association with identity and position


% Sequence & Periodicity
seq_autocorr

% Intrinsic timescales
acf_intrinsic_timescales
acf_intrinsic_timescales2