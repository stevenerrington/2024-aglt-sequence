% ========================================================================
% Script: bootstrap_acf_timescale.m
%
% Purpose:
%   Bootstrap spike-count autocorrelation functions (ACFs) from baseline
%   spike activity, estimate neuronal timescales (tau) by fitting an
%   exponential decay model, and compare timescales between auditory and
%   frontal neuron populations.
%
% High-Level Workflow:
%   1. Load / cache results:
%        - If cached ACFs (acf_boot_out.mat) do not exist, compute them.
%   2. Per neuron, per bootstrap iteration:
%        - Sample 25 valid trials (with replacement) meeting condition:
%          cond_label == 'nonviol' & rewardOnset_ms finite.
%        - Bin spikes in baseline window [-1000, 0] ms using 50 ms bins.
%        - Compute across-trial Pearson correlations for all bin pairs.
%        - Collapse correlation matrix by absolute lag -> 1D ACF vs lag (ms).
%        - Store ACF in acf_boot_out(neuron, lag, boot).
%   3. Save bootstrapped ACFs to disk (or load if already computed).
%   4. For each area (auditory, frontal):
%        - Select neurons_in = intersection of area idx and nonzero_neurons.
%        - For each bootstrap: average ACF across selected neurons (fit range
%          50–500 ms) and fit exponential decay A*exp(-lag/tau)+C via lsqcurvefit.
%        - Store tau in tau_test(:, area_i).
%   5. Visualize bootstrap tau distributions (histograms).
%   6. Compute frontal–auditory tau difference distribution, 95% CI, and a
%      simple two-tailed p-estimate based on sign of differences.
%   7. Plot tau difference using gramm (jitter + percentile summary).
%
% Dependencies:
%   MATLAB Optimization Toolbox (lsqcurvefit), gramm (for plotting), figuren()
%   utility (custom, assumed in path), Statistics Toolbox (randsample, histcounts).
%
% Author: Steven Errington   Date: 2025-07-21
% ========================================================================


%% Extract sound aligned Spike Density Function (SDF)
if ~exist(fullfile(dirs.root,'data','acf_boot_out.mat')) % Check if the output file already exists
    fprintf('Extracting bootstrapped autocorrelations \n')

    % Loop over all neurons in spike_log
    for neuron_i = 1:size(spike_log,1)
        fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

        % Clear variables to avoid contamination between neurons
        clear sdf_in event_table_in spikeTimesCell validTrials

        % ---------------- Load Data ----------------
        % Load spike data for this neuron (contains raster/trials/spike times)
        sdf_in = load(fullfile(dirs.root,'data','spike', ...
            [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

        % Load corresponding event table for this session
        event_table_in = load(fullfile(dirs.mat_data, ...
            [spike_log.session{neuron_i} '.mat']), 'event_table');

        % Extract cell array of spike times for each trial
        spikeTimesCell = sdf_in.raster.trialStart;

        % ---------------- Analysis Parameters ----------------
        baselineStart = -1000;  % Start of baseline window relative to event (ms)
        baselineEnd = 0;        % End of baseline window (ms)
        binSize = 50;           % Bin width (ms)
        edges = baselineStart:binSize:baselineEnd;  % Bin edges from -1000 to 0
        nBins = length(edges)-1;  % Number of bins (should be 20 for 1000 ms / 50 ms)
        delta = binSize;        % Lag step size (same as bin size, 50 ms)

        % Identify valid trials (condition = 'nonviol' and reward onset exists)
        validTrials_all = find(strcmp(event_table_in.event_table.cond_label, 'nonviol') & ...
            ~isnan(event_table_in.event_table.rewardOnset_ms));

        % ---------------- Bootstrap Loop ----------------
        for boot_i = 1:nboot

            try
                % Sample 25 trials with replacement from the valid trials
                validTrials = randsample(validTrials_all, 25, true);

                % ---------------- Compute Spike Counts ----------------
                nTrials = length(validTrials);
                spikeCountsMat = zeros(nTrials, nBins); % Matrix: trials x bins

                % Bin spike counts for baseline period for each trial
                for t = 1:nTrials
                    spikeTimes = spikeTimesCell{validTrials(t)};
                    spikeCountsMat(t,:) = histcounts(spikeTimes, edges);
                end

                % ---------------- Remove Bins with Zero Mean Firing ----------------
                meanFiringPerBin = mean(spikeCountsMat, 1); % Mean spike count per bin across trials
                nonzeroBins = meanFiringPerBin > 0;         % Logical mask for bins with nonzero mean

                if any(nonzeroBins == 0)
                    zero_fire_flag(neuron_i,boot_i) = 1; % Flag if any bins have zero firing
                else
                    zero_fire_flag(neuron_i,boot_i) = 0;
                end

                % ---------------- Compute Across-Trial Autocorrelation ----------------
                rhoMatrix = nan(nBins, nBins); % Matrix to store correlations between bins

                for idx_i = 1:nBins
                    for idx_j = 1:nBins
                        x = spikeCountsMat(:,idx_i); % Spike counts in bin i across trials
                        y = spikeCountsMat(:,idx_j); % Spike counts in bin j across trials

                        % Compute Pearson correlation only if both bins have variance
                        if std(x) > 0 && std(y) > 0
                            rhoMatrix(idx_i,idx_j) = corr(x, y); % Correlation across trials
                        end
                    end
                end

                % ---------------- Collapse by Lag ----------------
                % Compute average correlation for each lag (distance |i-j| between bins)
                maxLag = nBins-1;                % Maximum lag = number of bins - 1
                lags = (0:maxLag) * delta;       % Lags in ms
                acf = nan(1, maxLag+1);          % Autocorrelation function (ACF)

                for lag = 0:maxLag
                    vals = [];                   % Collect correlations for this lag
                    for idx_i = 1:(nBins-lag)
                        idx_j = idx_i + lag;
                        vals(end+1) = rhoMatrix(idx_i,idx_j); % Diagonal elements at lag
                    end
                    acf(lag+1) = mean(vals, 'omitnan'); % Average across pairs for this lag
                end

                % Store ACF for this neuron and bootstrap iteration
                acf_boot_out(neuron_i,:,boot_i) = acf;
            catch
                % If error occurs, fill with NaNs
                acf_boot_out(neuron_i,:,boot_i) = nan(1, length(lags));
            end
        end
    end

    % Save bootstrapped ACF results
    save(fullfile(dirs.root,'data','acf_boot_out.mat'),'acf_boot_out','-v7.3')
else
    fprintf('Loading bootstrapped acf data \n')
    load(fullfile(dirs.root,'data','acf_boot_out.mat'));
end

%% Select neurons and fit ACF decay

neurons_in = frontal_neuron_idx; % Default to frontal neurons

for area_i = 1:2
    switch area_i
        case 1
            area = 'auditory';
            neurons_in = []; neurons_in = intersect(auditory_neuron_idx,nonzero_neurons);
        case 2
            area = 'frontal';
            neurons_in = []; neurons_in = intersect(frontal_neuron_idx,nonzero_neurons);
    end

    for boot_i = 1:nboot
        % ---------------- Aggregate Across Neurons ----------------
        fitRange = lags >= 50 & lags <= 500;  % Use lags between 50 and 500 ms for fitting
        acf_input = nanmean(acf_boot_out(neurons_in,fitRange,boot_i)); % Average ACF for neuron class

        % ---------------- Fit Exponential Decay ----------------
        % Model: A * exp(-lag/tau) + C
        expDecay = @(b,x) b(1)*exp(-x/b(2)) + b(3);          % 3 parameters: amplitude, tau, offset
        expDecayOffset = @(b,x) b(1) * (exp(-x/b(2)) + b(3));% Alternative model form

        b0 = [0.5, 100, 0];  % Initial guess: amplitude=0.5, tau=100 ms, offset=0
        fitRange = lags >= 50 & lags <= 500;  % Fit only in this range

        % Perform nonlinear least squares fitting
        b_fit = lsqcurvefit(expDecay, b0, lags(fitRange), acf_input, [], []);
        decay_fitted = expDecay(b_fit, lags(fitRange));

        % Store tau (the timescale) for this bootstrap
        tau_test(boot_i,area_i) = b_fit(2);
    end
end

% Plot tau histograms for auditory (1) and frontal (2) neurons
figuren('Renderer', 'painters', 'Position', [500 300 500 300]);
histogram(tau_test(:,1),0:10:500,'LineStyle','none')
histogram(tau_test(:,2),0:10:500,'LineStyle','none')
xlabel([char([0xD835 0xDF0F]) ' (ms)']) % Math font tau
ylabel('N_{ bootstrap iterations}')

%% Diff data: Compare auditory vs frontal taus

tau_diff = tau_test(:,2)-tau_test(:,1);          % Difference distribution
ci = prctile(tau_diff, [2.5 97.5]);              % 95% CI for difference
p_val = 2 * min(mean(tau_diff >= 0), mean(tau_diff <= 0));  % Two-tailed p-value

% Plot difference distribution with jittered points and 95% percentile summary
clear fig_tau_diff
figure('Renderer', 'painters', 'Position', [760 350 200 200]);
fig_tau_diff=gramm('x',repmat({'Difference'},length(tau_diff),1),'y',tau_diff);
fig_tau_diff.geom_jitter('alpha',0.1);
fig_tau_diff.stat_summary('type','95percentile','geom',{'black_point','black_errorbar'});
fig_tau_diff.geom_hline(); % Zero line
fig_tau_diff.set_names('y','Difference (ms)');
fig_tau_diff.axe_property('YLim',[-100 300]);
fig_tau_diff.draw();
