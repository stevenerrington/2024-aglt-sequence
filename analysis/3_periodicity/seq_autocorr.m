rng(1,"twister"); % Set random number generator seed for reproducibility using 'twister' algorithm
acorr_time = 1000+[0:2665]; % Define time window for autocorrelation (ms)

%%
sound_times = [0, 563, 1126, 1689, 2252]; % Define sound onset times (ms) relative to some reference
acorr_neuron = [];
% Loop through each neuron
parfor neuron_i = 1:size(spike_log,1)
    % Display progress for the current neuron
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1));

    % Load the spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

    % Load the event table for the current session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

    % Extract spike density function (SDF) for 'nonviol' condition
    nonviol_sdf = [];
    nonviol_sdf = sdf_in.sdf.sequenceOnset(strcmp(event_table_in.event_table.cond_label, 'nonviol'), :);

    % Calculate baseline firing rate mean and standard deviation over a pre-stimulus period (-200 ms to 0 ms)
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, [1000+[-200:0]])));
    baseline_fr_std = nanstd(nanmean(nonviol_sdf(:, [1000+[-200:0]])));

    % Compute and smooth normalized SDF for the current neuron across time points
    seq_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);

    seq_sdf_out_temp = seq_sdf_out(neuron_i,:);
    [acorr_neuron(neuron_i,:), ~] = xcorr(seq_sdf_out_temp(:,acorr_time), 'coeff');

end

%% Run cluster average autocorrelation

% Compute population-averaged autocorrelation for auditory and frontal neurons
[acorr_aud_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.auditory.all,acorr_time)), 'coeff');
[acorr_frontal_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.frontal.all,acorr_time)), 'coeff');


figuren
plot(lags, acorr_aud_pos)
plot(lags, acorr_frontal_pos)
xlim([-1000 1000])


%%
maxLag = 2000; % in ms or samples
x = nanmean(seq_sdf_out(neuron_class.auditory.all,acorr_time));
y = nanmean(seq_sdf_out(neuron_class.frontal.all,acorr_time));

[C,lags] = xcorr(x-mean(x), y-mean(y), maxLag, 'coeff'); % normalized

numShuffles = 100;
% Optional: trial shuffle control
C_shuffle = nan(numShuffles, length(lags));
for s = 1:numShuffles
    y_shuff = y(randperm(length(y))); % shuffle trials
    C_shuffle(s,:) = xcorr(x-mean(x), y_shuff-mean(y_shuff), maxLag, 'coeff');
end
C_corrected = C - mean(C_shuffle,1); % shift predictor

figuren;
plot(lags,C)

CI = prctile(C_shuffle, [2.5 97.5], 1);
figure;
plot(lags,C_corrected); hold on;
plot(lags, CI(1,:), '--r');
plot(lags, CI(2,:), '--r');
xlabel('Lag (samples)'); ylabel('Cross-corr'); title('Auditory ↔ Frontal');
