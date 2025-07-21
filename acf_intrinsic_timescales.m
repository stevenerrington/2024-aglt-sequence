
% Remaining questions
% - smoothing

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
    validTrials = find(strcmp(event_table_in.event_table.cond_label, 'nonviol') & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms));

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
        zero_fire_flag(neuron_i,1) = 1;
    else
        zero_fire_flag(neuron_i,1) = 0;
    end

    % ---------------- Compute Across-Trial Autocorrelation ----------------
    rhoMatrix = nan(nBins, nBins); % Stores correlation between bins (i,j)
   
    for i = 1:nBins
        for j = 1:nBins
            x = spikeCountsMat(:,i); % Spike counts in bin i across trials
            y = spikeCountsMat(:,j); % Spike counts in bin j across trials
            
            % Compute Pearson correlation only if both bins have variance
            if std(x) > 0 && std(y) > 0
                rhoMatrix(i,j) = corr(x, y); % Correlation across trials
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
        for i = 1:(nBins-lag)
            j = i + lag;
            vals(end+1) = rhoMatrix(i,j); % Diagonal elements at lag
        end
        acf(lag+1) = mean(vals, 'omitnan'); % Average across pairs for this lag
    end

    % Store ACF for this neuron
    acf_out(neuron_i,:) = acf;
    rhoMatrix_out{neuron_i} = rhoMatrix;
end


nonzero_neurons = find(zero_fire_flag == 0);

%%
figure('Renderer', 'painters', 'Position', [100 100 600 600]);
    clear plot_intrinsic_timescale

for area_i = 1:2
    switch area_i
        case 1
            area = 'auditory';
            neurons_in = []; neurons_in = intersect(auditory_neuron_idx,nonzero_neurons);
        case 2
            area = 'frontal';
            neurons_in = []; neurons_in = intersect(frontal_neuron_idx,nonzero_neurons);
    end

    % ---------------- Aggregate Across Neurons ----------------
    fitRange = lags >= 50 & lags <= 500;  % Use lags between 50 and 500 ms for fitting
    acf_input = nanmean(acf_out(neurons_in,fitRange)); % Average ACF for neuron class

    % ---------------- Fit Exponential Decay ----------------
    % Model: A * exp(-lag/tau) + C
    expDecay = @(b,x) b(1)*exp(-x/b(2)) + b(3);
    expDecayOffset = @(b,x) b(1) * (exp(-x/b(2)) + b(3));

    b0 = [0.5, 100, 0];  % Initial guess: amplitude=0.5, tau=100 ms, offset=0
    fitRange = lags >= 50 & lags <= 500;  % Fit only in this range

    % Perform nonlinear least squares fitting
    b_fit = lsqcurvefit(expDecay, b0, lags(fitRange), acf_input, [], []);
    decay_fitted = expDecay(b_fit, lags(fitRange));


    % -----------------
    % Plot intrinsic timescale
    plot_intrinsic_timescale(area_i,1) = gramm('x', lags(fitRange), 'y', acf_out(neurons_in,fitRange));
    plot_intrinsic_timescale(area_i,1).stat_summary('geom',{'point','line','errorbar'});
    plot_intrinsic_timescale.axe_property('YLim',[0 0.1]);

    plot_intrinsic_timescale(area_i,2) = gramm('x', lags(fitRange), 'y', decay_fitted);
    plot_intrinsic_timescale(area_i,2).geom_line('alpha',1);
    plot_intrinsic_timescale.axe_property('YLim',[0 0.1]);

    tau_intrinsic_timescale_fit(area_i) = b_fit(2);

    SStot = sum((acf_input-mean(acf_input)).^2);                            % Total Sum-Of-Squares
    SSres = sum((acf_input(:)-decay_fitted(:)).^2);                         % Residual Sum-Of-Squares
    r2_intrinsic_timescale_fit(area_i) = 1-SSres/SStot;                                    % R^2
end


plot_intrinsic_timescale.draw;
snapnow;




%%
figure('Renderer', 'painters', 'Position', [100 100 300 650]);

for area_i = 1:2
    switch area_i
        case 1
            area = 'auditory';
            neurons_in = []; neurons_in = intersect(auditory_neuron_idx,nonzero_neurons);
        case 2
            area = 'frontal';
            neurons_in = []; neurons_in = intersect(frontal_neuron_idx,nonzero_neurons);
    end

    rhoMatrix_average = [];
    % Plot autocorrelation heatmap
    for neuron_i = 1:length(neurons_in)
        rhoMatrix_average(:,:,neuron_i) = rhoMatrix_out{neurons_in(neuron_i)};
    end

    plot_rhoMatrix_average = [];
    plot_rhoMatrix_average = nanmean(rhoMatrix_average,3);
    plot_rhoMatrix_average(plot_rhoMatrix_average > 0.99) = NaN;

    subplot(2,1,area_i)
    h(area_i) = heatmap(plot_rhoMatrix_average);
    h(area_i).Colormap = hot;
    h(area_i).ColorLimits = [0 0.08];
    h(area_i).MissingDataColor = [78 122 199]./255;      % optional: show NaNs as white
    h(area_i).MissingDataLabel = 'NaN';        % optional label in colorbar
    h(area_i).GridVisible = 'off';

end




