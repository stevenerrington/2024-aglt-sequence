function get_laminar_profile(electrode_ref, lfp, agl_t_event)

electrode_list = electrode_ref;
clear signal_out

%% Restructure data for spectrolaminar analysis

fprintf('Restructuring LFP data into the appropriate format   \n')
% - for each electrode
for electrode_i = 1:length(electrode_list)
    electrode_idx = electrode_list(electrode_i);
    n_trials = size(lfp.(['lfp_' int2str(electrode_i)]),1);

    % - across each trial
    for trial_i = 1:n_trials

        % - get the LFP data for the given
        signal_in = lfp.(['lfp_' int2str(electrode_idx)])(trial_i,:);

        % - save the relevant data in an output array for future use
        signal_out(electrode_i,:,trial_i) = signal_in; % nchans x trialtime x ntrials


    end
end

%% Find valid trials
valid_trial_idx = ~isnan(agl_t_event.stimulusOnset_ms);

%% Get signal power (normalized and non-normalized; required FieldTrip)
fprintf('Calculating signal power   \n')

clear signal_normalized signal_normalized
[signal_normalized, signal_nonnormalized] = lfp_to_normpower(signal_out(:,:,valid_trial_idx));

%% Run the main spectrolaminar extraction analysis code
fprintf('Running FLIP toolbox   \n')

% FLIP (Spectrolaminar) toolbox
[startinglowfreq,endinglowfreq,...
    startinghighfreq,endinghighfreq,...
    goodnessvalue,superficialchannel,deepchannel,...
    highfreqmaxchannel,lowfreqmaxchannel,...
    crossoverchannel] = ...
    FLIPAnalysis(signal_normalized,0:size(signal_normalized,1)-1,1:size(signal_normalized,2),1);

% Cross-contact correlation
fprintf('Running xcorr toolbox   \n')
xcontact_corr = D_LFPCORR_BASIC(signal_out, 512, true);

%% Plot figure
fprintf('Generating figure   \n')

clear color_heatmap
color_heatmap = cbrewer2('seq', 'GnBu', 100);

figuren('Renderer', 'painters', 'Position', [100 100 1200 300]); hold on;
subplot(1,5,[1 2]); hold on;
imagesc(signal_normalized);set(gca, 'YDir', 'reverse');
xlim([1 size(signal_normalized,2)]); ylim([1 size(signal_normalized,1)]);
xlabel('Frequency (Hz)');ylabel('Channel Number');
colormap(color_heatmap)
hline(crossoverchannel,'k')
cb=colorbar; ylabel(cb, 'Relative Power');
% clim([0.5 1.0]); 

subplot(1,5,3); hold on;
plot(mean(signal_normalized(:,10:19),2), 1:length(electrode_list), 'b', 'LineWidth', 2); % alpha/beta
plot(mean(signal_normalized(:,75:150),2), 1:length(electrode_list), 'r', 'LineWidth', 2); % gamma
set(gca, 'ydir', 'reverse')
xlim([0 1]); ylim([1 length(electrode_list)]); xlabel('Relative Power'); ylabel('Channel Number');
hline(crossoverchannel,'k')

subplot(1,5,[4 5]); hold on;
imagesc(xcontact_corr);set(gca, 'YDir', 'reverse');
xlim([1 size(xcontact_corr,2)]); ylim([1 size(xcontact_corr,1)]);
hline(crossoverchannel,'k')
vline(crossoverchannel,'k')
cb=colorbar; ylabel(cb, 'Correlation');
%clim([0.8 1.0]);

end
