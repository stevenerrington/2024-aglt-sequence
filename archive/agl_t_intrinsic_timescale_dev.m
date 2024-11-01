acorr_time = 1000+[0:2665];
spike_train = nanmean(seq_sdf_out(neuron_class.auditory.all,acorr_time))


% Step 1: Define parameters
dt = 1;  % time step in ms
maxLag = 500;  % max lag in ms to compute autocorrelation (adjust as needed)

% Step 2: Compute autocorrelation of spike train
[acf, lags] = xcorr(spike_train, maxLag, 'coeff'); % normalized ACF

% Step 3: Fit exponential decay to the positive lags of ACF
positiveLags = lags(lags >= 49); 
acf_positive = acf(lags >= 49);
acf_positive = acf_positive - min(acf_positive);

% Define exponential decay function for fitting
exp_decay = @(b, x) b(1) * exp(-x / b(2));  % b(1) is the amplitude, b(2) is tau

% Initial guess for [amplitude, tau]
initial_guess = [1, 500];  % adjust initial guess as needed

% Perform curve fitting
fit_params = nlinfit(positiveLags, acf_positive, exp_decay, initial_guess);

% Extract intrinsic timescale (tau)
tau = fit_params(2);  % tau in ms

% Plot the autocorrelation and the fitted curve
figure;
plot(positiveLags, acf_positive, 'b', 'LineWidth', 1.5); hold on;
plot(positiveLags, exp_decay(fit_params, positiveLags), 'r--', 'LineWidth', 1.5);
xlabel('Lag (ms)');
ylabel('Autocorrelation');
legend('Autocorrelation', 'Exponential Fit');
title(['Intrinsic Timescale (\tau) = ', num2str(tau), ' ms']);
xlim([0 maxLag + 50])

