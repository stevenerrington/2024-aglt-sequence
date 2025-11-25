%% ================================================================
% Intrinsic Timescale Analysis
% Computes spike-count autocorrelation (ACF) across trials for each neuron,
% fits an exponential decay to estimate intrinsic timescale, and
% performs bootstrapped analyses for population-level statistics.
% ================================================================

%% ------------------ Parameters ------------------
baselineStart = -1000;          % Baseline window start (ms)
baselineEnd   = 0;              % Baseline window end (ms)
binSize       = 50;             % Bin size (ms)
edges         = baselineStart:binSize:baselineEnd;
nBins         = length(edges) - 1;
delta         = binSize;        % Lag step (ms)

lags    = (0:nBins-1) * delta;  % Lag vector
fitMask = lags >= 50 & lags <= 500;  % Lags used for exponential fit

%% ------------------ Preallocate outputs ------------------
nNeurons        = size(spike_log,1);
acf_out         = nan(nNeurons, nBins);    % Autocorrelation per neuron
rhoMatrix_out   = cell(nNeurons,1);       % Spike-count correlation matrices
zero_fire_flag  = zeros(nNeurons,1);      % Flag neurons with zero firing in any bin

%% ================================================================
% Loop Over Neurons: Compute Spike-Count ACF
% ================================================================
parfor neuron_i = 1:nNeurons
    fprintf('Processing neuron %i of %i\n', neuron_i, nNeurons);

    %% ---------- Load Data ----------
    sdf_in = load(fullfile(dirs.root,'data','spike', ...
        sprintf('%s_%s.mat', spike_log.session{neuron_i}, spike_log.unitDSP{neuron_i})));
    evt = load(fullfile(dirs.mat_data, ...
        sprintf('%s.mat', spike_log.session{neuron_i})), 'event_table');

    spikeTimesCell = sdf_in.raster.trialStart;

    %% ---------- Select Valid Trials ----------
    validTrials = find(strcmp(evt.event_table.cond_label, 'nonviol') & ...
                       ~isnan(evt.event_table.rewardOnset_ms));
    nTrials = length(validTrials);

    %% ---------- Bin Spike Counts ----------
    spikeCounts = zeros(nTrials, nBins);
    for t = 1:nTrials
        spikeCounts(t,:) = histcounts(spikeTimesCell{validTrials(t)}, edges);
    end

    %% ---------- Drop Bins with Zero Mean Firing ----------
    meanFiring = mean(spikeCounts,1);
    nonzeroBins = meanFiring > 0;
    zero_fire_flag(neuron_i) = any(~nonzeroBins);

    %% ---------- Compute Spike-Count Correlation Matrix ----------
    rhoMat = nan(nBins, nBins);
    for i = 1:nBins
        xi = spikeCounts(:,i);
        if std(xi)==0, continue; end
        for j = 1:nBins
            yj = spikeCounts(:,j);
            if std(yj)==0, continue; end
            rhoMat(i,j) = corr(xi, yj);
        end
    end
    rhoMatrix_out{neuron_i} = rhoMat;

    %% ---------- Collapse Correlation by Lag ----------
    acf = nan(1,nBins);
    for lag = 0:nBins-1
        diagVals = diag(rhoMat, lag);
        acf(lag+1) = mean(diagVals, 'omitnan');
    end
    acf_out(neuron_i,:) = acf;
end

nonzero_neurons = find(~zero_fire_flag);

%% ================================================================
% Plot: Intrinsic Timescale Fits
% ================================================================
figure('Renderer','painters','Position',[100 100 600 600]);

areas   = {'auditory','frontal'};
idxList = {auditory_neuron_idx, frontal_neuron_idx};

for area_i = 1:2
    neurons_in = intersect(idxList{area_i}, nonzero_neurons);

    %% Aggregate across neurons
    acf_mean = nanmean(acf_out(neurons_in, fitMask));

    %% Fit exponential decay: A*exp(-lag/tau) + C
    expDecay = @(b,x) b(1)*exp(-x/b(2)) + b(3);  % Model: amplitude, tau, offset
    b0 = [0.5, 100, 0];                          % Initial guess

    b_fit = lsqcurvefit(expDecay, b0, lags(fitMask), acf_mean);
    decay_fitted = expDecay(b_fit, lags(fitMask));

    %% Plot ACF and fitted decay using gramm
    plot_intrinsic_timescale(area_i,1) = gramm('x', lags(fitMask), 'y', acf_out(neurons_in,fitMask));
    plot_intrinsic_timescale(area_i,1).stat_summary('geom',{'point','line','errorbar'});
    plot_intrinsic_timescale(area_i,1).axe_property('YLim',[0 0.1]);

    plot_intrinsic_timescale(area_i,2) = gramm('x', lags(fitMask), 'y', decay_fitted);
    plot_intrinsic_timescale(area_i,2).geom_line();
    plot_intrinsic_timescale(area_i,2).axe_property('YLim',[0 0.1]);

    tau_intrinsic_timescale_fit(area_i) = b_fit(2);

    %% Compute R²
    SS_tot = sum((acf_mean - mean(acf_mean)).^2);
    SS_res = sum((acf_mean - decay_fitted).^2);
    r2_intrinsic_timescale_fit(area_i) = 1 - SS_res/SS_tot;
end

plot_intrinsic_timescale.draw;
snapnow;

%% ================================================================
% Plot: Average Spike-Count Correlation Matrices
% ================================================================
figure('Renderer','painters','Position',[100 100 300 650]);
for area_i = 1:2
    neurons_in = intersect(idxList{area_i}, nonzero_neurons);

    rhoAvg = nanmean(cat(3, rhoMatrix_out{neurons_in}), 3);
    rhoAvg(rhoAvg > 0.99) = NaN;

    subplot(2,1,area_i)
    h(area_i) = heatmap(rhoAvg);
    h(area_i).Colormap = hot;
    h(area_i).ColorLimits = [0 0.08];
    h(area_i).MissingDataColor = [78 122 199]/255;
    h(area_i).MissingDataLabel = 'NaN';
    h(area_i).GridVisible = 'off';
end

%% ================================================================
% Bootstrapped Spike-Count Autocorrelation (ACF)
% ================================================================
rng(1, "twister");                  % Reproducibility
acorr_time = 1000 + (0:2665);       % Time points for ACF calculation

sound_times = [0, 563, 1126, 1689, 2252]; % Sound onset times

acorr_neuron = [];

parfor neuron_i = 1:nNeurons
    fprintf('Neuron %i of %i \n', neuron_i, nNeurons);

    % Load neuron data
    sdf_in = load(fullfile(dirs.root,'data','spike', ...
        [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
    event_table_in = load(fullfile(dirs.mat_data, ...
        [spike_log.session{neuron_i} '.mat']), 'event_table');

    % Extract SDF for 'nonviol' trials
    nonviol_sdf = sdf_in.sdf.sequenceOnset(strcmp(event_table_in.event_table.cond_label, 'nonviol'), :);

    % Baseline normalization (-200 to 0 ms)
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:,1000+(-200:0))));
    baseline_fr_std  = nanstd(nanmean(nonviol_sdf(:,1000+(-200:0))));

    seq_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);

    % Compute autocorrelation
    seq_sdf_out_temp = seq_sdf_out(neuron_i,:);
    [acorr_neuron(neuron_i,:), ~] = xcorr(seq_sdf_out_temp(:,acorr_time), 'coeff');
end

%% ================================================================
% Population-Averaged Autocorrelation
% ================================================================
[acorr_aud_pos, lags]    = xcorr(nanmean(seq_sdf_out(auditory_neuron_idx,acorr_time)), 'coeff');
[acorr_frontal_pos, lags] = xcorr(nanmean(seq_sdf_out(frontal_neuron_idx,acorr_time)), 'coeff');

figure;
plot(lags, acorr_aud_pos);
hold on;
plot(lags, acorr_frontal_pos);
xlim([-1000 1000]);

%% ================================================================
% Bootstrapped Peak Autocorrelation
% ================================================================
nboot   = 1000;
nSample = 500;
peak_autocorr = nan(nboot,2);

for boot_i = 1:nboot
    acorr_aud_pos_boot    = xcorr(nanmean(seq_sdf_out(randsample(auditory_neuron_idx, nSample, true), acorr_time)), 'coeff');
    acorr_frontal_pos_boot = xcorr(nanmean(seq_sdf_out(randsample(frontal_neuron_idx, nSample, true), acorr_time)), 'coeff');

    % Find first positive peak
    try
        [~,peak_idx_aud] = findpeaks(acorr_aud_pos_boot,'MinPeakProminence',0.20);
        peak_autocorr(boot_i,1) = lags(peak_idx_aud(find(lags(peak_idx_aud)>0,1)));
    catch
        peak_autocorr(boot_i,1) = NaN;
    end
    try
        [~,peak_idx_frontal] = findpeaks(acorr_frontal_pos_boot,'MinPeakProminence',0.20);
        peak_autocorr(boot_i,2) = lags(peak_idx_frontal(find(lags(peak_idx_frontal)>0,1)));
    catch
        peak_autocorr(boot_i,2) = NaN;
    end
end

% Plot histogram of peaks
figure;
histogram(peak_autocorr(:,1),0:10:750,'LineStyle','None'); hold on;
histogram(peak_autocorr(:,2),0:10:750,'LineStyle','None');

fprintf('Auditory valid: %.2f%%\n', 100*(1 - mean(isnan(peak_autocorr(:,1)))));
fprintf('Frontal valid: %.2f%%\n', 100*(1 - mean(isnan(peak_autocorr(:,2)))));

%% ================================================================
% Bootstrapped Exponential Fitting of ACF
% ================================================================
if ~exist(fullfile(dirs.root,'data','acf_boot_out.mat'),'file')
    fprintf('Extracting bootstrapped autocorrelations \n');

    for neuron_i = 1:nNeurons
        fprintf('Neuron %i of %i \n', neuron_i, nNeurons);

        sdf_in = load(fullfile(dirs.root,'data','spike', ...
            [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
        event_table_in = load(fullfile(dirs.mat_data, ...
            [spike_log.session{neuron_i} '.mat']), 'event_table');

        spikeTimesCell = sdf_in.raster.trialStart;
        validTrials_all = find(strcmp(event_table_in.event_table.cond_label, 'nonviol') & ...
            ~isnan(event_table_in.event_table.rewardOnset_ms));

        for boot_i = 1:nboot
            try
                validTrials = randsample(validTrials_all, length(validTrials_all), true);
                nTrials = length(validTrials);
                spikeCountsMat = zeros(nTrials, nBins);

                for t = 1:nTrials
                    spikeCountsMat(t,:) = histcounts(spikeTimesCell{validTrials(t)}, edges);
                end

                % Remove bins with zero mean firing
                meanFiringPerBin = mean(spikeCountsMat,1);
                zero_fire_flag(neuron_i,boot_i) = any(meanFiringPerBin==0);

                % Compute correlation matrix
                rhoMatrix = nan(nBins, nBins);
                for iBin = 1:nBins
                    for jBin = 1:nBins
                        x = spikeCountsMat(:,iBin);
                        y = spikeCountsMat(:,jBin);
                        if std(x)>0 && std(y)>0
                            rhoMatrix(iBin,jBin) = corr(x,y);
                        end
                    end
                end

                % Collapse by lag
                acf = nan(1,nBins);
                for lag = 0:nBins-1
                    acf(lag+1) = mean(diag(rhoMatrix, lag), 'omitnan');
                end

                acf_boot_out(neuron_i,:,boot_i) = acf;
            catch
                acf_boot_out(neuron_i,:,boot_i) = nan(1,nBins);
            end
        end
    end
    save(fullfile(dirs.root,'data','acf_boot_out.mat'),'acf_boot_out','-v7.3')
else
    fprintf('Loading bootstrapped ACF data \n');
    load(fullfile(dirs.root,'data','acf_boot_out.mat'));
end

%% ================================================================
% Fit exponential decay to bootstrapped ACFs
% ================================================================

tau_test = nan(nboot,2);
for area_i = 1:2
    switch area_i
        case 1
            neurons_in = intersect(auditory_neuron_idx, nonzero_neurons);
        case 2
            neurons_in = intersect(frontal_neuron_idx, nonzero_neurons);
    end

    for boot_i = 1:nboot


        % lags from baseline binning
        lags_bins = (0:nBins-1) * delta;

        % Only use lags that exist in your ACF
        fitMask_bins = lags_bins >= 50 & lags_bins <= 500;

        acf_input = nanmean(acf_boot_out(neurons_in, fitMask_bins, boot_i));

        b_fit = lsqcurvefit(expDecay, b0, lags_bins(fitMask_bins), acf_input, [], []);

        expDecay = @(b,x) b(1)*exp(-x/b(2)) + b(3);
        b0 = [0.5, 100, 0];

        b_fit = lsqcurvefit(expDecay, b0, lags_bins(fitMask_bins), acf_input, [], []);
        tau_test(boot_i,area_i) = b_fit(2);
    end
end

% Plot histogram of bootstrapped taus
figure('Renderer','painters','Position',[500 300 500 300]);
histogram(tau_test(:,1),0:10:500,'LineStyle','none'); hold on;
histogram(tau_test(:,2),0:10:500,'LineStyle','none');
xlabel('\tau (ms)');
ylabel('Bootstrap iterations');

% Compute 95% CI
ci_auditory = prctile(tau_test(:,1), [2.5 50 97.5]);
ci_frontal  = prctile(tau_test(:,2), [2.5 50 97.5]);

% Difference distribution
tau_diff = tau_test(:,2)-tau_test(:,1);
ci_diff = prctile(tau_diff,[2.5 50 97.5]);
p_val = 2*min(mean(tau_diff>=0), mean(tau_diff<=0));

% Plot difference with jittered points
figure('Renderer','painters','Position',[760 350 200 200]);
fig_tau_diff = gramm('x', repmat({'Difference'}, length(tau_diff),1), 'y', tau_diff);
fig_tau_diff.geom_jitter('alpha',0.1);
fig_tau_diff.stat_summary('type','95percentile','geom',{'black_point','black_errorbar'});
fig_tau_diff.geom_hline();
fig_tau_diff.set_names('y','Difference (ms)');
fig_tau_diff.axe_property('YLim',[-100 300]);
fig_tau_diff.draw();
