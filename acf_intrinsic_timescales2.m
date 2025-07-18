
%% Extract sound aligned Spike Density Function (SDF)
if ~exist(fullfile(dirs.root,'data','acf_boot_out.mat'))
    fprintf('Extracting bootstrapped autocorrelations \n')

    % Loop over all neurons in spike_log
    for neuron_i = 1:size(spike_log,1)
        fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

        % Clear variables to avoid contamination between neurons
        clear sdf_in event_table_in spikeTimesCell validTrials

        % ---------------- Load Data ----------------
        % Load spike data for this neuron
        sdf_in = load(fullfile(dirs.root,'data','spike', ...
            [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

        % Load corresponding event table for session
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

        for boot_i = 1:nboot

            try
                validTrials = randsample(validTrials_all, 25, true);
                % ---------------- Compute Spike Counts ----------------
                nTrials = length(validTrials);
                spikeCountsMat = zeros(nTrials, nBins); % Matrix: trials x bins

                % Bin spike counts for baseline period
                for t = 1:nTrials
                    spikeTimes = spikeTimesCell{validTrials(t)};
                    spikeCountsMat(t,:) = histcounts(spikeTimes, edges);
                end

                % ---------------- Remove Bins with Zero Mean Firing ----------------
                meanFiringPerBin = mean(spikeCountsMat, 1); % Mean spike count per bin across trials
                nonzeroBins = meanFiringPerBin > 0;        % Logical mask for bins with nonzero mean

                if any(nonzeroBins == 0)
                    zero_fire_flag(neuron_i,boot_i) = 1;
                else
                    zero_fire_flag(neuron_i,boot_i) = 0;
                end

                % ---------------- Compute Across-Trial Autocorrelation ----------------
                rhoMatrix = nan(nBins, nBins); % Stores correlation between bins (i,j)

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
                % Compute average correlation for each lag (|i-j| bins)
                maxLag = nBins-1;                % Maximum lag = number of bins - 1
                lags = (0:maxLag) * delta;       % Lags in ms
                acf = nan(1, maxLag+1);          % Autocorrelation function

                for lag = 0:maxLag
                    vals = [];                   % Collect correlations for this lag
                    for idx_i = 1:(nBins-lag)
                        idx_j = idx_i + lag;
                        vals(end+1) = rhoMatrix(idx_i,idx_j); % Diagonal elements at lag
                    end
                    acf(lag+1) = mean(vals, 'omitnan'); % Average across pairs for this lag
                end

                % Store ACF for this neuron
                acf_boot_out(neuron_i,:,boot_i) = acf;
            catch
                acf_boot_out(neuron_i,:,boot_i) = nan(1, length(lags));
            end
        end
    end

    save(fullfile(dirs.root,'data','acf_boot_out.mat'),'acf_boot_out','-v7.3')
else
    fprintf('Loading bootstrapped acf data \n')
    load(fullfile(dirs.root,'data','acf_boot_out.mat'));
end

%%

neurons_in = frontal_neuron_idx;

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
        expDecayOld = @(b,x) b(1)*exp(-x/b(2)) + b(3);
        expDecayOffset = @(b,x) b(1) * (exp(-x/b(2)) + b(3));

        b0 = [0.5, 100, 0];  % Initial guess: amplitude=0.5, tau=100 ms, offset=0
        fitRange = lags >= 50 & lags <= 500;  % Fit only in this range

        % Perform nonlinear least squares fitting
        b_fit = lsqcurvefit(expDecayOffset, b0, lags(fitRange), acf_input, [], []);
        decay_fitted = expDecayOffset(b_fit, lags(fitRange));

        tau_test(boot_i,area_i) = b_fit(2);
    end
end

figuren;
histogram(tau_test(:,1),0:10:500,'LineStyle','none')
histogram(tau_test(:,2),0:10:500,'LineStyle','none')


figuren;
histogram(tau_test(:,2)-tau_test(:,1),-500:10:500)
