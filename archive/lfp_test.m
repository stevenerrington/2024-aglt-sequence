session_idx = 10;
session = ephysLog.session{session_idx};
% Load the corresponding session data from the specified directory
session_data = load(fullfile(dirs.mat_data, session));

% Define frequencies of interest for analysis
freqs_in = [4, 8; % Theta
    10, 15;   % Alpha
    15, 30;   % Beta
    40, 70];  % Gamma

filters = design_filters(freqs_in);

timewin = -1000:5000;
lfp_data_in = []; lfp_data_in = session_data.lfp;
trials_in = []; trials_in = find(~isnan(session_data.event_table.stimulusOnset_ms) & ...
    ~isnan(session_data.event_table.rewardOnset_ms) & strcmp(session_data.event_table.cond_label,'nonviol'));

alignEvent = 'stimulusOnset_ms';
alignTimes = session_data.event_table.(alignEvent) + session_audio_latency{session_idx};
[phase_analysis_out] = get_lfp_measures(lfp_data_in, alignTimes(trials_in), filters, timewin);

plot_freq = 'gamma';
plot_electrode = 26;
plot_xlim = [-1000 4000];
plot_vline = [0 563 1126 1689 2252];

figure('Renderer', 'painters', 'Position', [100 100 500 800]); hold on;
subplot(6,1,1); box off
plot(timewin, nanmean(phase_analysis_out.([plot_freq '_lfp']){plot_electrode}))
xlim(plot_xlim); vline(plot_vline, 'k')
title('Raw Frequency LFP')

subplot(6,1,2); box off
plot(timewin, nanmean(phase_analysis_out.([plot_freq '_amplitude']){plot_electrode}))
xlim(plot_xlim); vline(plot_vline, 'k')
title('Amplitude')

subplot(6,1,3); box off
plot(timewin, phase_analysis_out.([plot_freq '_plv']){plot_electrode})
xlim(plot_xlim); vline(plot_vline, 'k')
title('Phase-locking value')

subplot(6,1,[4 5 6]); box off
imagesc('XData',timewin, 'YData', 1:size(phase_analysis_out.([plot_freq '_phase']){plot_electrode},1),...
    'CData', phase_analysis_out.([plot_freq '_phase']){plot_electrode})
xlim(plot_xlim); ylim([1 size(phase_analysis_out.([plot_freq '_phase']){plot_electrode},1)]); vline(plot_vline, 'k')
xlabel(['Time from ' alignEvent])
title('Phase')

suptitle([session ' - ' plot_freq])

