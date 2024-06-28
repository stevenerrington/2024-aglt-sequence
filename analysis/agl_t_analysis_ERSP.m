%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- AGLt ERSP analysis -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%} 

% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find event information to align data on
aligntime = event_table.stimulusOnset_ms;

clear trials*
trials_viol = find(~isnan(aligntime) & strcmp(event_table.cond_label,'viol'));
trials_nonviol = find(~isnan(aligntime) & strcmp(event_table.cond_label,'nonviol'));
trials_input = find(~isnan(aligntime));

% Run alignment algorithms
ops.timewin = -1000:5000;
ops.freq = [2 200];
lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

% Restructure electrode data into a ch x time x trial format
electrode_list = [1:32];
% - for each electrode
for electrode_i = 1:length(electrode_list)
    electrode_idx = electrode_list(electrode_i);
    n_trials = size(lfp_aligned.(['lfp_' int2str(electrode_i)]),1);

    % - across each trial
    for trial_i = 1:n_trials

        % - get the LFP data for the given
        signal_in = lfp_aligned.(['lfp_' int2str(electrode_idx)])(trial_i,:);

        % - save the relevant data in an output array for future use
        signal_out(electrode_i,:,trial_i) = signal_in; % nchans x trialtime x ntrials


    end
end

data = signal_out(:,:,trials_input);

% EEGlab analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters and configuration -----------------------------------------
srate =  1000; % Sampling rate (Hz)
epochmin = ops.timewin(1)/1000; % Epoch start (sec)
epochmax = ops.timewin(end)/1000; % Epoch end (sec)
basemin = -500; % Baseline window start (ms)
basemax = 0; % Baseline window end (ms)

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
clear data EEG

% - Structure data for EEGlab -------------------------------------------
% (adapted from YK earlier matlab code)
data = signal_out(:,:,trials_input);
EEG = pop_importdata('dataformat', 'array', 'data', 'data', 'srate',srate, 'nbchan',16);
EEG = eeg_checkset(EEG);

% - Define epochs (although data is already aligned)
EEG = pop_importepoch(EEG, in, { 'Epoch', 'stim'}, 'latencyfields',{ 'stim'}, 'timeunit',1, 'headerlines',0);
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  'stim'  }, [epochmin         epochmax], 'newname', 'Level epochs', 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );

% - Perform baseline correction
EEG = pop_rmbase( EEG, [basemin    0]);
EEG = eeg_checkset( EEG );

% Run EEGlab time-frequency analyses -------------------------------------
chan_i = 3;

clear ersp itc
fprintf('Running analysis on channel %i \n', chan_i)
[ersp,itc,powbase,times,freqs,erspboot,itcboot,alltfX] = pop_newtimef(EEG, ...
    1, chan_i, [EEG.xmin EEG.xmax]*srate, [3 0.7], 'maxfreq',maxfreq, 'freqs',freq_range,'padratio', padratio, ...
    'plotphase', 'off', 'timesout', outtimes, 'alpha', alpha_val, 'naccu', 200, 'baseboot',1,'rmerp','off', ...
    'erspmax', maxersp, 'plotersp','off', 'plotitc','off','baseline',[basemin basemax],'marktimes',0);



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

