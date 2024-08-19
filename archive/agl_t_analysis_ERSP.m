%%
%{
///////////////////////////////////////////////////////////////////////////
----- AGLt ERSP analysis -----------------------------------------
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

%% Analysis ///////////////////////////////////////////////////////////////
tic
for session_i = 1:size(ephysLog,1)
    load(fullfile(dirs.mat_data,ephysLog.session{session_i}))

    % Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find event information to align data on
    aligntime = event_table.stimulusOnset_ms;

    clear trials*
    trials_viol = find(~isnan(aligntime) & strcmp(event_table.cond_label,'viol'));
    trials_nonviol = find(~isnan(aligntime) & strcmp(event_table.cond_label,'nonviol'));
    trials_input = find(~isnan(aligntime));

    % Run alignment algorithms
    ops.timewin = -1000:5000;
    ops.freq = [2 200];
    [~, signal_out] = get_lfp_aligned(lfp,aligntime,ops);

    data = signal_out(:,:,trials_nonviol);

    % EEGlab analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters and configuration -----------------------------------------
    srate =  1000; % Sampling rate (Hz)
    epochmin = ops.timewin(1)/1000; % Epoch start (sec)
    epochmax = ops.timewin(end)/1000; % Epoch end (sec)
    basemin = -500; % Baseline window start (ms)
    basemax = 0; % Baseline window end (ms)

    clear nTr in
    nTr=size(data,3); % Number of trials
    in(:,1)=[1:nTr]'; % Dummy variable for trial n
    in(:,2)=abs(epochmin)*ones(nTr,1); % Dummy variable for trial n

    freq_range=[2.5 100]; % Frequency range for ERSP analysis
    maxfreq = max(freq_range); % Max frequency
    padratio = 2; % Pad ratio
    outtimes = 1500; % Out times
    alpha_val = 0.01; % Alpha
    maxersp = 6;


    % Setup data ------------------------------------------------------------
    clear EEG

    % - Structure data for EEGlab -------------------------------------------
    % (adapted from YK earlier matlab code)
    EEG = pop_importdata('dataformat', 'array', 'data', 'data', 'srate',srate, 'nbchan',32);
    EEG = eeg_checkset(EEG);

    % - Define epochs (although data is already aligned)
    EEG = pop_importepoch(EEG, in, { 'Epoch', 'stim'}, 'latencyfields',{ 'stim'}, 'timeunit',1, 'headerlines',0);
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'stim'  }, [epochmin         epochmax], 'newname', 'Level epochs', 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );

    % - Perform baseline correction
    EEG = pop_rmbase( EEG, [basemin    0]);
    EEG = eeg_checkset( EEG );

    clear ersp_out
    % Run EEGlab time-frequency analyses -------------------------------------
    for chan_i = 1:32

        clear ersp itc powbase times freqs erspboot itcboot alltfX

        fprintf('Running analysis on channel %i \n', chan_i)
        [ersp,itc,powbase,times,freqs,erspboot,itcboot,alltfX] = pop_newtimef(EEG, ...
            1, chan_i, [EEG.xmin EEG.xmax]*srate, [3 0.7], 'maxfreq',maxfreq, 'freqs',freq_range,'padratio', padratio, ...
            'plotphase', 'off', 'alpha', alpha_val, 'naccu', 200, 'baseboot',1,'rmerp','off', ...
            'erspmax', maxersp, 'plotersp','off', 'plotitc','off','baseline',[basemin basemax],'marktimes',0);

        ersp_out{chan_i}.ersp = ersp;
        ersp_out{chan_i}.powbase = powbase;
        ersp_out{chan_i}.times = times;
        ersp_out{chan_i}.freqs = freqs;
        ersp_out{chan_i}.erspboot = erspboot;
        ersp_out{chan_i}.itcboot = itcboot;
        ersp_out{chan_i}.alltfX = alltfX;

    end
        
    
    save(fullfile('C:\KIKUCHI-LOCAL\data\ephys\ersp',[ephysLog.session{session_i} '-ersp.mat']),'ersp_out','-v7.3')

end

toc



% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colorscale = 'parula';

figure;
subplot(2,1,1)
imagesc('XData',ops.timewin,'YData',freqs,'CData',ersp)
xlim([-500 4000]); ylim([min(freqs) 100])
vline(0, 'k-');
colorbar; colormap(colorscale)

subplot(2,1,2)
imagesc('XData',ops.timewin,'YData',freqs,'CData',abs(itc))
xlim([-500 4000]); ylim([min(freqs) 100]); clim([0 0.5])
vline(0, 'k-');
colorbar;

