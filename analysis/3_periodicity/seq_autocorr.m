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
[acorr_aud_pos, lags] = xcorr(nanmean(seq_sdf_out(auditory_neuron_idx,acorr_time)), 'coeff');
[acorr_frontal_pos, lags] = xcorr(nanmean(seq_sdf_out(frontal_neuron_idx,acorr_time)), 'coeff');


figuren;
plot(lags, acorr_aud_pos)
plot(lags, acorr_frontal_pos)
xlim([-1000 1000])

%%

nboot = 1000;
nSample = 500;


for boot_i = 1:nboot

    [acorr_aud_pos_boot] = xcorr(nanmean(seq_sdf_out(randsample(auditory_neuron_idx, nSample, true),acorr_time)), 'coeff');
    [acorr_frontal_pos_boot] = xcorr(nanmean(seq_sdf_out(randsample(frontal_neuron_idx, nSample, true),acorr_time)), 'coeff');

    clear peak_idx_aud peak_idx_frontal
    [~,peak_idx_aud] = findpeaks(acorr_aud_pos_boot,'MinPeakProminence',0.20);
    [~,peak_idx_frontal] = findpeaks(acorr_frontal_pos_boot,'MinPeakProminence',0.20);

    try
        peak_autocorr(boot_i,1) = lags(peak_idx_aud(find(lags(peak_idx_aud) > 0,1)));
    catch
        peak_autocorr(boot_i,1) = NaN;
    end

    try
        peak_autocorr(boot_i,2) = lags(peak_idx_frontal(find(lags(peak_idx_frontal) > 0,1)));
    catch
        peak_autocorr(boot_i,2) = NaN;
    end


end


figuren;
histogram(peak_autocorr(:,1),0:10:750,'LineStyle','None')
histogram(peak_autocorr(:,2),0:10:750,'LineStyle','None')

disp(100-(sum(isnan(peak_autocorr(:,1)))./length(peak_autocorr))*100)
disp(100-(sum(isnan(peak_autocorr(:,2)))./length(peak_autocorr))*100)


mode(peak_autocorr(:,1))
mode(peak_autocorr(:,2))

mode(peak_autocorr(peak_autocorr(:,1) < 500,1))
mode(peak_autocorr(peak_autocorr(:,2) < 500,1))


%%
