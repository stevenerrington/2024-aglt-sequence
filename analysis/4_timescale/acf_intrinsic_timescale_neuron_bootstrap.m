%% Fig 4 intrinsic-timescale: a genuine neuron-resampling bootstrap CI
%
% acf_intrinsic_timescales.m's bootstrap section resamples TRIALS within
% each neuron (with replacement, 1000 times) to build acf_boot_out, but
% then always averages that trial-resampled ACF across the SAME, FIXED,
% full set of neurons (neurons_in) on every single bootstrap iteration --
% which neurons contribute to the population mean never changes. That
% cancels out within-neuron trial-sampling noise but never lets
% between-neuron variability (some neurons have a longer/shorter timescale
% than others) propagate into the confidence interval on tau, or on the
% frontal-vs-auditory difference in tau. A CI that only reflects trial
% noise, not neuron-to-neuron variability, will generally be too narrow --
% worth re-doing properly given how close the reported difference already
% is to the significance threshold (73 ms, CI [6, 150], p = 0.034).
%
% Fix: resample WHICH NEURONS enter the population-average ACF on each
% bootstrap iteration (with replacement, drawing as many neurons as are
% actually available in that area -- standard nonparametric bootstrap
% practice, matching every other neuron-resampling bootstrap already used
% elsewhere in this analysis, e.g. pca_seq_main.m, pca_ramp_structure_
% test.m), refitting the exponential decay fresh on each draw. This is
% scoped ONLY to the missing-neuron-resampling issue: the per-neuron ACF
% computation (spike-count correlation matrix, collapsed by lag), the
% zero-firing neuron-exclusion criterion, the exponential decay model,
% and the (unconstrained) lsqcurvefit call are all kept identical to
% acf_intrinsic_timescales.m so this is a like-for-like comparison. The
% other two issues flagged for this analysis (no positivity constraint on
% tau; zero-firing exclusion drops a whole neuron rather than one bad
% bin) are deliberately left untouched here and can be addressed
% separately.
%
% Because resampling neurons only requires each neuron's own REAL
% (non-trial-resampled) ACF, this needs just one pass over the raw spike
% data (not 1000 passes) -- the per-neuron ACF is computed once, then
% resampled-and-refit 1000 times at negligible extra cost.
%
% Assumes the following are already in the workspace: spike_log, dirs,
% auditory_neuron_idx, frontal_neuron_idx.

clear acf_out zero_fire_flag

rng(1,'twister')

%% ---- Parameters (identical to acf_intrinsic_timescales.m) ----
baselineStart = -1000;
baselineEnd   = 0;
binSize       = 50;
edges         = baselineStart:binSize:baselineEnd;
nBins         = length(edges) - 1;
delta         = binSize;

lags    = (0:nBins-1) * delta;
fitMask = lags >= 50 & lags <= 500;

expDecay = @(b,x) b(1)*exp(-x/b(2)) + b(3);
b0       = [0.5, 100, 0];

n_boot = 1000;

%% ---- Step 1: real (non-bootstrapped) per-neuron ACF -- one pass over the data ----

nNeurons       = size(spike_log,1);
acf_out        = nan(nNeurons, nBins);
zero_fire_flag = zeros(nNeurons,1);

for neuron_i = 1:nNeurons

    if mod(neuron_i,100) == 0
        fprintf('Neuron %i of %i\n', neuron_i, nNeurons);
    end

    try
        sdf_in = load(fullfile(dirs.root, 'data', 'spike', ...
            sprintf('%s_%s.mat', spike_log.session{neuron_i}, spike_log.unitDSP{neuron_i})));
        evt = load(fullfile(dirs.mat_data, ...
            sprintf('%s.mat', spike_log.session{neuron_i})), 'event_table');
    catch
        zero_fire_flag(neuron_i) = 1;
        continue
    end

    spikeTimesCell = sdf_in.raster.trialStart;
    validTrials = find(strcmp(evt.event_table.cond_label, 'nonviol') & ...
                       ~isnan(evt.event_table.rewardOnset_ms));
    nTrials = length(validTrials);

    if nTrials < 2
        zero_fire_flag(neuron_i) = 1;
        continue
    end

    spikeCounts = zeros(nTrials, nBins);
    for t = 1:nTrials
        spikeCounts(t,:) = histcounts(spikeTimesCell{validTrials(t)}, edges);
    end

    meanFiring  = mean(spikeCounts,1);
    nonzeroBins = meanFiring > 0;
    zero_fire_flag(neuron_i) = any(~nonzeroBins);

    rhoMat = nan(nBins, nBins);
    for i = 1:nBins
        xi = spikeCounts(:,i);
        if std(xi) == 0, continue; end
        for j = 1:nBins
            yj = spikeCounts(:,j);
            if std(yj) == 0, continue; end
            rhoMat(i,j) = corr(xi, yj);
        end
    end

    acf = nan(1,nBins);
    for lag = 0:nBins-1
        acf(lag+1) = mean(diag(rhoMat, lag), 'omitnan');
    end
    acf_out(neuron_i,:) = acf;
end

nonzero_neurons = find(~zero_fire_flag);
fprintf('Neurons with valid (non-zero-firing) ACF: %d / %d\n', numel(nonzero_neurons), nNeurons);

%% ---- Step 2: point-estimate tau per area (matches acf_intrinsic_timescales.m's first fit) ----

decoding_problems = {'Auditory', auditory_neuron_idx; 'Frontal', frontal_neuron_idx};

tau_point  = nan(1,2);
neurons_in_area = cell(1,2);

for area_i = 1:2
    neurons_in = intersect(decoding_problems{area_i,2}, nonzero_neurons);
    neurons_in_area{area_i} = neurons_in;

    acf_mean = nanmean(acf_out(neurons_in, fitMask));
    b_fit = lsqcurvefit(expDecay, b0, lags(fitMask), acf_mean, [], []);
    tau_point(area_i) = b_fit(2);

    fprintf('%s: n = %d neurons, point-estimate tau = %.1f ms\n', ...
        decoding_problems{area_i,1}, numel(neurons_in), tau_point(area_i));
end

%% ---- Step 3: neuron-resampling bootstrap ----

tau_boot = nan(n_boot,2);

for area_i = 1:2

    label      = decoding_problems{area_i,1};
    neurons_in = neurons_in_area{area_i};
    nSample    = numel(neurons_in);   % resample the actual sample size -- standard nonparametric bootstrap

    fprintf('\n=== %s: neuron-resampling bootstrap (n_boot = %d, nSample = %d) ===\n', label, n_boot, nSample);

    for boot_i = 1:n_boot
        boot_idx = randsample(neurons_in, nSample, true);
        acf_boot_mean = nanmean(acf_out(boot_idx, fitMask));

        try
            b_fit = lsqcurvefit(expDecay, b0, lags(fitMask), acf_boot_mean, [], []);
            tau_boot(boot_i,area_i) = b_fit(2);
        catch
            tau_boot(boot_i,area_i) = nan;
        end
    end

    ci = prctile(tau_boot(:,area_i), [2.5 50 97.5]);
    fprintf('tau: median = %.1f ms, 95%% CI [%.1f, %.1f]\n', ci(2), ci(1), ci(3));
end

%% ---- Step 4: frontal vs auditory, properly neuron-resampled ----

fprintf('\nFrontal vs auditory tau (neuron-resampling bootstrap):\n');
bootstrap_compare(tau_boot(:,2), tau_boot(:,1));

tau_diff = tau_boot(:,2) - tau_boot(:,1);
ci_diff = prctile(tau_diff, [2.5 50 97.5]);
fprintf('Difference (frontal - auditory): median = %.1f ms, 95%% CI [%.1f, %.1f]\n', ...
    ci_diff(2), ci_diff(1), ci_diff(3));
fprintf('(For reference, the manuscript-reported difference was 73 ms, 95%% CI [6, 150], p = 0.034,\n');
fprintf(' from a bootstrap that never resampled neurons.)\n');

%% ---- Step 5: plot ----

figure('Renderer','painters','Position',[100 100 600 400]);

subplot(1,2,1); hold on
area_colors = [0.2 0.4 0.7; 0.8 0.3 0.2];
areas = {'Auditory','Frontal'};
for area_i = 1:2
    histogram(tau_boot(:,area_i), 0:10:500, 'FaceColor', area_colors(area_i,:), ...
        'FaceAlpha', 0.5, 'EdgeColor','none');
end
xlabel('\tau (ms)'); ylabel('Bootstrap iterations (neuron-resampled)');
legend(areas); box off
title('Neuron-resampled \tau distributions')

subplot(1,2,2); hold on
bar(1, median(tau_diff,'omitnan'), 0.5, 'FaceColor',[0.5 0.5 0.5]);
errorbar(1, median(tau_diff,'omitnan'), median(tau_diff,'omitnan')-ci_diff(1), ...
    ci_diff(3)-median(tau_diff,'omitnan'), 'k', 'LineWidth',1.2);
yline(0,'k--');
set(gca,'XTick',1,'XTickLabel',{'Frontal - Auditory'});
ylabel('\Delta\tau (ms)'); box off
title('Difference, neuron-resampled CI')
