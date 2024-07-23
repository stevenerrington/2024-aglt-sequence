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

% Set extraction directories
dirs.py_eventtable = 'C:\KIKUCHI-LOCAL\data\ephys\py_conv\event_table';
dirs.py_lfp = 'C:\KIKUCHI-LOCAL\data\ephys\py_conv\lfp';
dirs.spike_times = 'C:\KIKUCHI-LOCAL\data\ephys\py_conv\spike_times';
dirs.spike_waves = 'C:\KIKUCHI-LOCAL\data\ephys\py_conv\spike_waves';

% For each session, loop through, extract the data and save it
% independently
for session_i = 1:size(ephysLog,1)
    fprintf('Converting session %i of %i to text/csv | %s \n', session_i, size(ephysLog,1), ephysLog.session{session_i})
    load(fullfile(dirs.mat_data,ephysLog.session{session_i}))
    
    % Export events to csv
    writetable(event_table,fullfile(dirs.py_eventtable,[ephysLog.session{session_i} '-events.csv']),...
        'WriteVariableNames',true)

    session_neurons_dsp = spike_log.unitDSP(strcmp(spike_log.session,ephysLog.session{session_i}));
    session_neurons_wav = spike_log.unitWAV(strcmp(spike_log.session,ephysLog.session{session_i}));

    % Export spike times and waveforms to txt
    for neuron_i = 1:length(session_neurons_dsp)
        writematrix(spikes.time.(session_neurons_dsp{neuron_i})',fullfile(dirs.spike_times,[ephysLog.session{session_i} '-' session_neurons_dsp{neuron_i} '-spk.txt']))
        writematrix(spikes.waveform.(session_neurons_wav{neuron_i}),fullfile(dirs.spike_waves,[ephysLog.session{session_i} '-' session_neurons_dsp{neuron_i} '-wav.txt']))
    end

    % Export spike times and waveforms to txt
    writematrix(lfp ,fullfile(dirs.py_lfp,[ephysLog.session{session_i} '-lfp.txt']))
end


