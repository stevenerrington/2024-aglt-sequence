% Define frequencies of interest for analysis
freqs_in = [4, 8; % Theta
    10, 15;   % Alpha
    15, 30;   % Beta
    40, 70];  % Gamma
filters = design_filters(freqs_in);

timewin = -1000:5000;

% Retrieve the session name from ephysLog_hpc based on the specified area and session index
for session_idx = 1:size(ephysLog, 1)
    session = ephysLog.session{session_idx};
    fprintf('Session %i of %i \n', session_idx, size(ephysLog,1));
    % Load the corresponding session data from the specified directory
    session_data = load(fullfile(dirs.mat_data, session));

    clear phase_analysis_out alignTimes lfp_data_in trials_in

    lfp_data_in = []; lfp_data_in = session_data.lfp;

    trials_in = []; trials_in = find(~isnan(session_data.event_table.rewardOnset_ms) & strcmp(session_data.event_table.cond_label,'nonviol'));

    alignEvent = 'rewardOnset_ms';
    alignTimes = session_data.event_table.(alignEvent) + session_audio_latency{session_idx};

    [phase_analysis_out] = get_lfp_measures(lfp_data_in, alignTimes(trials_in), filters, timewin);

    save(fullfile(dirs.root,'data','phase_analysis',['phase_analysis_out_session' int2str(session_idx) '.mat']),'phase_analysis_out','-v7.3')
end


%%
session_idx = 60;
load(fullfile(dirs.root,'data','phase_analysis',['phase_analysis_out_session' int2str(session_idx) '.mat']))

elec_ch = 28;

figuren
imagesc(timewin, 1:size(phase_analysis_out.theta_phase{elec_ch},1), phase_analysis_out.theta_phase{elec_ch})
colorscale_phase = cbrewer('div', 'RdBu', 11);  % Green color scale for backward transitions
colormap(colorscale_phase); colorbar
xlim([-1000 1000]); vline(0,'k-')
ylim([1 size(phase_analysis_out.theta_phase{elec_ch},1)])




















%%
fs = 1000; % Hz, adjust to your sampling rate
timewin = -1000:5000;  % in ms
baseline_idx = find((timewin >= -1000 & timewin < 0) | (timewin >= 3500 & timewin < 5000));
sequence_idx = find(timewin >= 0 & timewin < 2500);

parfor session_idx = 1:size(ephysLog, 1)

    session = ephysLog.session{session_idx};
    fprintf('Session %i of %i \n', session_idx, size(ephysLog,1));

    % Load the corresponding session data from the specified directory
    session_data = load(fullfile(dirs.mat_data, session));
    lfp_data_in = []; lfp_data_in = session_data.lfp;
    trials_in = []; trials_in = find(~isnan(session_data.event_table.stimulusOnset_ms) & ~isnan(session_data.event_table.rewardOnset_ms) & strcmp(session_data.event_table.cond_label,'nonviol'));

    alignEvent = 'stimulusOnset_ms';
    alignTimes = round(session_data.event_table.(alignEvent) + session_audio_latency{session_idx});
    alignTimes = alignTimes;

    [~, lfp_trials] = get_lfp_aligned(session_data.lfp,alignTimes,ops);

    for ch = 1:32

        % grab data (concatenate trials along time)
        baseline = reshape(lfp_trials(ch, baseline_idx, trials_in), 1, []);
        sequence = reshape(lfp_trials(ch, sequence_idx, trials_in), 1, []);

        baseline = baseline(~isnan(baseline));
        sequence = sequence(~isnan(sequence));
        
        % Welch parameters
        win = hanning(2*fs); % 2 s window
        noverlap = length(win)/2;
        nfft = [];

        [Pf_base{session_idx}(ch,:),f] = pwelch(baseline, win, noverlap, nfft, fs);
        [Pf_seq{session_idx}(ch,:), ~] = pwelch(sequence, win, noverlap, nfft, fs);

    end
end

%%

parfor session_idx = 1:size(ephysLog, 1)
    fprintf('Session %i of %i \n', session_idx, size(ephysLog,1));
    probe1_area = ephysLog.area_label_sec{session_idx};
    frontal_power = []; auditory_power = []; frontal_plv = []; auditory_plv = [];


    switch probe1_area
        case 'auditory'
            for electrode_i = 1:16
                auditory_pf_base(:,:,session_idx) = Pf_base{session_idx}(1:16,:);
                auditory_pf_seq(:,:,session_idx) = Pf_seq{session_idx}(1:16,:);
                frontal_pf_base(:,:,session_idx) = Pf_base{session_idx}(17:32,:);
                frontal_pf_seq(:,:,session_idx) = Pf_seq{session_idx}(17:32,:);
            end
        case 'frontal'
            for electrode_i = 1:16
                auditory_pf_base(:,:,session_idx) = Pf_base{session_idx}(17:32,:);
                auditory_pf_seq(:,:,session_idx) = Pf_seq{session_idx}(17:32,:);
                frontal_pf_base(:,:,session_idx) = Pf_base{session_idx}(1:16,:);
                frontal_pf_seq(:,:,session_idx) = Pf_seq{session_idx}(1:16,:);
            end
    end
end


%%

figuren;
subplot(1,2,1); hold on
plot(f, nanmean(nanmean(auditory_pf_base,3)))
plot(f, nanmean(nanmean(auditory_pf_seq,3)))
xlim([0 40]); title('Auditory')
subplot(1,2,2); hold on
plot(f, nanmean(nanmean(frontal_pf_base,3)))
plot(f, nanmean(nanmean(frontal_pf_seq,3)))
xlim([0 40]); title('Frontal')
legend({'baseline','sequence'})


%%

theta_idx = find(f > freqs_in(1,1) & f < freqs_in(1,2));
alpha_idx = find(f > freqs_in(2,1) & f < freqs_in(2,2));
beta_idx = find(f > freqs_in(3,1) & f < freqs_in(3,2));
gamma_idx = find(f > freqs_in(4,1) & f < freqs_in(4,2));


clear auditory_pwr_change frontal_pwr_change 
auditory_pwr_change = {[], [], [], []};
frontal_pwr_change = {[], [], [], []};
for session_idx = 1:size(ephysLog, 1)
    for freq_i = 1:4
        auditory_pwr_change{freq_i} = [auditory_pwr_change{freq_i}; 100*(nanmean(auditory_pf_seq(:, f > freqs_in(freq_i,1) & f < freqs_in(freq_i,2), session_idx),2)./nanmean(auditory_pf_base(:, f > freqs_in(freq_i,1) & f < freqs_in(freq_i,2), session_idx),2)-1)];
        frontal_pwr_change{freq_i} = [frontal_pwr_change{freq_i}; 100*(nanmean(frontal_pf_seq(:, f > freqs_in(freq_i,1) & f < freqs_in(freq_i,2), session_idx),2)./nanmean(frontal_pf_base(:, f > freqs_in(freq_i,1) & f < freqs_in(freq_i,2), session_idx),2)-1)];
    end
end


power_labels = repmat([repmat({'1_Theta'}, 1264, 1); repmat({'2_Alpha'}, 1264, 1); repmat({'3_Beta'}, 1264, 1); repmat({'4_Gamma'}, 1264, 1)],2,1);
power_data = [auditory_pwr_change{1}; auditory_pwr_change{2}; auditory_pwr_change{3}; auditory_pwr_change{4};...
    frontal_pwr_change{1}; frontal_pwr_change{2}; frontal_pwr_change{3}; frontal_pwr_change{4}];
power_area = [repmat({'Auditory'}, 1264*4, 1); repmat({'Frontal'}, 1264*4, 1)];

% Plot bootstrap classification accuracy
figure('Renderer', 'painters', 'Position', [100 100 400 300]);
clear power_figure
power_figure(1,1) = gramm('x', power_area, 'y', power_data, 'color', power_labels);
power_figure(1,1).stat_summary('dodge',0.7,'geom',{'bar','black_errorbar'});
power_figure(1,1).axe_property('YLim',[-25 125]);
power_figure(1,1).set_names('X', 'Area'); % Set the x-axis label name
power_figure(1,1).set_names('Y', '% change from baseline'); % Set the y-axis label name
power_figure.draw;
