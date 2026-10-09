%% Per-neuron intrinsic timescale (tau), not just the population-averaged fit
%
% Every timescale analysis so far (acf_intrinsic_timescales.m, and the
% neuron-resampling bootstrap fix) fits ONE exponential decay to the ACF
% averaged across all neurons in an area, then bootstraps that population-
% level fit. That's the right way to get a robust area-level tau and CI
% (and it's the number that held up: 71.3 ms frontal-vs-auditory
% difference, CI [7.5, 138.5], p = 0.026) -- but it never produces a tau
% value for any individual neuron, which is needed for anything at the
% single-unit level (distributions, correlating tau with other per-neuron
% properties, identifying which specific neurons carry the effect, etc.).
%
% This script fits the same exponential decay model, over the same lag
% window (50-500 ms) and the same per-neuron ACF computation as the
% neuron-resampling bootstrap script, but to EACH NEURON'S OWN ACF
% individually rather than to an area-pooled average. Single-neuron ACFs
% are necessarily much noisier than a several-hundred-neuron average (they
% come from correlating spike counts across that one neuron's own trials,
% typically a few dozen), so an unconstrained fit can return degenerate
% values (runaway or negative tau). Two things are added here that were
% deliberately NOT added to the population-level fit: loose but real
% bounds on the fit (tau in [1, 3000] ms; amplitude/offset in [-5, 5], far
% wider than the ACF values themselves ever go, just to stop the
% optimizer diverging) and a per-neuron R^2 (goodness of the exponential
% fit to that neuron's own ACF), so obviously bad fits can be filtered
% rather than silently included. Neurons whose fitted tau sits exactly on
% the upper bound are flagged separately -- that's not a real 3000 ms
% timescale, it means the fit couldn't find a decay at all.
%
% If acf_out and nonzero_neurons are already in the workspace (e.g. from
% having just run acf_intrinsic_timescale_neuron_bootstrap.m), this reuses
% them directly instead of recomputing every neuron's ACF from raw spike
% data again. Otherwise it computes them fresh, identically.
%
% Inclusion criteria: matched to the standard already used for the same
% kind of analysis (autocorrelation exponential decay, Wasmuht/Zeisler
% approach) in the frontal ramping project -- at least 250 spikes (summed
% across the same trials/bins that go into that neuron's ACF), R^2 > 0.5,
% 0 < tau < 1000 ms, and a positive fitted amplitude (A > 0; a negative
% amplitude means the fit found a rise, not a decay, which isn't a
% meaningful "intrinsic timescale"). The fit itself still uses the loose
% [1, 3000] ms bound on tau (not [0, 1000] directly) so the optimizer
% isn't artificially prevented from finding a true fit outside the
% acceptance window -- the 1000 ms ceiling is then applied afterward, as
% an inclusion criterion, not a fitting constraint.
%
% Assumes the following are already in the workspace: spike_log, dirs,
% auditory_neuron_idx, frontal_neuron_idx.

rng(1,'twister')

%% ---- Parameters (identical to acf_intrinsic_timescales.m / the neuron-bootstrap fix) ----
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
lb       = [-5,    1, -5];    % loose bounds: stop the optimizer diverging, not a strong prior
ub       = [ 5, 3000,  5];

fit_opts = optimoptions('lsqcurvefit','Display','off');

%% ---- Step 1: per-neuron ACF (reuse if already computed) ----

if exist('acf_out','var') && exist('nonzero_neurons','var') && exist('neuron_total_spikes','var') ...
        && size(acf_out,1) == size(spike_log,1)

    fprintf('Reusing acf_out / nonzero_neurons / neuron_total_spikes already in the workspace.\n');

else

    fprintf('Computing per-neuron ACF and total spike count from raw spike data.\n');

    nNeurons          = size(spike_log,1);
    acf_out           = nan(nNeurons, nBins);
    zero_fire_flag    = zeros(nNeurons,1);
    neuron_total_spikes = nan(nNeurons,1);

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

        neuron_total_spikes(neuron_i) = sum(spikeCounts(:));   % total spikes across the same trials/bins feeding the ACF

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
end

fprintf('Neurons with valid ACF available for per-neuron fitting: %d / %d\n', ...
    numel(nonzero_neurons), size(spike_log,1));

%% ---- Step 2: fit tau to each neuron's own ACF individually ----

nNeurons = size(spike_log,1);
neuron_tau    = nan(nNeurons,1);
neuron_amp    = nan(nNeurons,1);
neuron_offset = nan(nNeurons,1);
neuron_r2     = nan(nNeurons,1);
neuron_at_bound = false(nNeurons,1);

for neuron_i = nonzero_neurons'

    y = acf_out(neuron_i, fitMask);
    x = lags(fitMask);

    if any(isnan(y))
        continue
    end

    try
        b_fit = lsqcurvefit(expDecay, b0, x, y, lb, ub, fit_opts);
    catch
        continue
    end

    yhat = expDecay(b_fit, x);
    ss_res = sum((y - yhat).^2);
    ss_tot = sum((y - mean(y)).^2);

    neuron_tau(neuron_i)    = b_fit(2);
    neuron_amp(neuron_i)    = b_fit(1);
    neuron_offset(neuron_i) = b_fit(3);
    neuron_r2(neuron_i)     = 1 - ss_res/ss_tot;
    neuron_at_bound(neuron_i) = abs(b_fit(2) - ub(2)) < 1;
end

%% ---- Step 3: apply the standard inclusion criteria (>=250 spikes, R^2>0.5, 0<tau<1000ms, A>0) ----

decoding_problems = {'Auditory', auditory_neuron_idx; 'Frontal', frontal_neuron_idx};

min_spikes  = 250;
r2_cut      = 0.5;
tau_lo      = 0;
tau_hi      = 1000;

pass_spikes = neuron_total_spikes >= min_spikes;
pass_r2     = neuron_r2 > r2_cut;
pass_tau    = neuron_tau > tau_lo & neuron_tau < tau_hi;
pass_amp    = neuron_amp > 0;
included    = pass_spikes & pass_r2 & pass_tau & pass_amp;

neuron_tau_table = table((1:nNeurons)', neuron_total_spikes, neuron_tau, neuron_amp, neuron_offset, neuron_r2, ...
    neuron_at_bound, included, ...
    'VariableNames', {'neuron_idx','total_spikes','tau','amplitude','offset','r2','tau_at_upper_bound','included'});

fprintf('\n=== Inclusion funnel (>=%d spikes, R^2>%.1f, %d<tau<%d ms, A>0) ===\n', min_spikes, r2_cut, tau_lo, tau_hi);
for area_i = 1:2
    label      = decoding_problems{area_i,1};
    neurons_in = intersect(decoding_problems{area_i,2}, nonzero_neurons);

    fprintf('\n%s (n = %d neurons with a valid ACF):\n', label, numel(neurons_in));
    fprintf('  pass spike-count (>=%d): %d\n', min_spikes, sum(pass_spikes(neurons_in)));
    fprintf('  pass R^2 (>%.1f):        %d\n', r2_cut, sum(pass_r2(neurons_in)));
    fprintf('  pass tau range (%d-%d ms): %d\n', tau_lo, tau_hi, sum(pass_tau(neurons_in)));
    fprintf('  pass amplitude (>0):     %d\n', sum(pass_amp(neurons_in)));

    sel = neurons_in(included(neurons_in));
    if isempty(sel)
        fprintf('  ALL CRITERIA: 0 neurons pass\n');
        continue
    end
    tau_sel = neuron_tau(sel);
    fprintf('  ALL CRITERIA: n = %d (%.0f%% of valid) | median tau = %.1f ms | IQR [%.1f, %.1f]\n', ...
        numel(sel), 100*numel(sel)/numel(neurons_in), median(tau_sel), prctile(tau_sel,25), prctile(tau_sel,75));
end

%% ---- Step 4: frontal vs auditory on the per-neuron distribution (all criteria applied) ----

aud_sel = intersect(auditory_neuron_idx, nonzero_neurons);
aud_sel = aud_sel(included(aud_sel));
fro_sel = intersect(frontal_neuron_idx, nonzero_neurons);
fro_sel = fro_sel(included(fro_sel));

fprintf('\nPer-neuron tau (all inclusion criteria applied), frontal vs auditory:\n');
fprintf('Wilcoxon rank-sum test (distribution-level, not the population-averaged-ACF fit):\n');
[p_ranksum,~,stats] = ranksum(neuron_tau(fro_sel), neuron_tau(aud_sel));
fprintf('  Frontal median = %.1f ms (n=%d) | Auditory median = %.1f ms (n=%d) | rank-sum p = %.4f\n', ...
    median(neuron_tau(fro_sel)), numel(fro_sel), median(neuron_tau(aud_sel)), numel(aud_sel), p_ranksum);

%% ---- Step 5: plot per-neuron tau distributions ----

figure('Renderer','painters','Position',[100 100 900 400]);

area_colors = [0.2 0.4 0.7; 0.8 0.3 0.2];
areas = {'Auditory','Frontal'};
sel_by_area = {aud_sel, fro_sel};

subplot(1,3,1); hold on
for area_i = 1:2
    histogram(neuron_tau(sel_by_area{area_i}), 0:20:1000, 'FaceColor', area_colors(area_i,:), ...
        'FaceAlpha', 0.5, 'EdgeColor','none', 'Normalization','probability');
end
xlabel('\tau (ms), per neuron'); ylabel('Proportion of neurons')
legend(areas); box off
title(sprintf('Per-neuron \\tau distribution (>=%d spikes, R^2>%.1f, %d<\\tau<%d, A>0)', min_spikes, r2_cut, tau_lo, tau_hi))

subplot(1,3,2); hold on
for area_i = 1:2
    scatter(neuron_r2(sel_by_area{area_i}), neuron_tau(sel_by_area{area_i}), 10, ...
        area_colors(area_i,:), 'filled', 'MarkerFaceAlpha', 0.4);
end
xlabel('Per-neuron fit R^2'); ylabel('\tau (ms)')
box off
title('Fit quality vs. tau (noisier fits at low R^2)')

subplot(1,3,3); hold on
for area_i = 1:2
    v = neuron_tau(sel_by_area{area_i});
    jitter_x = area_i + 0.15*(rand(size(v))-0.5);
    scatter(jitter_x, v, 8, area_colors(area_i,:), 'filled', 'MarkerFaceAlpha', 0.3);
    plot(area_i + [-0.2 0.2], median(v)*[1 1], 'k-', 'LineWidth', 2);
end
set(gca,'XTick',[1 2],'XTickLabel',areas)
ylabel('\tau (ms), per neuron'); box off
title('Median (black line) per area')

fprintf('\nPer-neuron results stored in neuron_tau_table (also as neuron_tau, neuron_r2, neuron_amp,\n');
fprintf('neuron_offset, neuron_total_spikes, neuron_at_bound, included vectors, one row per neuron in spike_log).\n');
