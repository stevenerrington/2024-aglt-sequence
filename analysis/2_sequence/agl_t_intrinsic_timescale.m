clear seq_sdf_out seq_sdf_out_rwd
% Loop through each neuron
parfor neuron_i = 1:size(spike_log,1)
    % Display progress for the current neuron
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1));

    % Load the spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

    % Load the event table for the current session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

    % Extract spike density function (SDF) for 'nonviolent' condition
    nonviol_sdf = [];
    nonviol_sdf = sdf_in.sdf.trialStart(strcmp(event_table_in.event_table.cond_label, 'nonviol'), :);

    % Calculate baseline firing rate mean and standard deviation over a pre-stimulus period (-200 ms to 0 ms)
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, [1000+[-200:0]])));
    baseline_fr_std = nanstd(nanmean(nonviol_sdf(:, [1000+[-200:0]])));

    % Compute and smooth normalized SDF for the current neuron across time points
    seq_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf(:,[1:1000])) - baseline_fr_mean) ./ baseline_fr_std, 50);
end

clear spike_train
spike_train = nanmean(seq_sdf_out(neuron_class.frontal.all,:));

bin_size = 50;  % Number of bins to average over
num_bins = floor(length(spike_train) / bin_size);  % Total groups

% Reshape and compute mean
spike_train_avg = mean(reshape(spike_train(1:num_bins * bin_size), bin_size, []), 1);

% Step 1: Define parameters
maxLag = 500;  % max lag in ms to compute autocorrelation (adjust as needed)

% Step 2: Compute autocorrelation of spike train
[acf, lags] = xcorr(spike_train, maxLag, 'coeff'); % normalized ACF

positiveLags = []; acf_positive = []; exp_decay = []; initial_guess = []; fit_params = []; tau = [];

% Step 3: Fit exponential decay to the positive lags of ACF
positiveLags = lags(lags >= 49);
acf_positive = acf(lags >= 49);
acf_positive = acf_positive - min(acf_positive);

% Define exponential decay function for fitting
exp_decay = @(b, x) b(1) * exp(-x / b(2));  % b(1) is the amplitude, b(2) is tau

% Initial guess for [amplitude, tau]
initial_guess = [1, 200];  % adjust initial guess as needed

% Perform curve fitting
fit_params = nlinfit(positiveLags, acf_positive, exp_decay, initial_guess);

% Extract intrinsic timescale (tau)
tau = fit_params(2);  % tau in ms



figuren;
subplot(2,1,1)
plot(-999:0, spike_train)
subplot(2,1,2); hold on
plot(positiveLags, exp_decay(fit_params, positiveLags), 'r-', 'LineWidth', 1.5);
plot(positiveLags, acf_positive, 'k-', 'LineWidth', 1.5);