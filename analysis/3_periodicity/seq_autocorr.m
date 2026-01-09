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


%% Plot sequence examples
figuren('Renderer','painters','Position',[100 100 300 500]);
subplot(4,1,1)
plot(-1000:5000, seq_sdf_out(1359,:))
xlim([-200 2750]); ylim([-2 8]); box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,1,2)
plot(-1000:5000, seq_sdf_out(1280,:))
xlim([-200 2750]); ylim([-2 20]); box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,1,3)
plot(-1000:5000, seq_sdf_out(721,:))
xlim([-200 2750]); ylim([-2 1]); box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,1,4)
plot(-1000:5000, seq_sdf_out(1352,:))
xlim([-200 2750]); ylim([-1 4]); box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')


%% Plot sequence cluster averages
figuren('Renderer','painters','Position',[100 100 600 500]);
% Auditory
subplot(4,2,1)
plot(-1000:5000, nanmean(seq_sdf_out(cluster_neurons.auditory{1},:)))
xlim([-200 2750]);  box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,2,3)
plot(-1000:5000, nanmean(seq_sdf_out(cluster_neurons.auditory{2},:)))
xlim([-200 2750]);  box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,2,5)
plot(-1000:5000, nanmean(seq_sdf_out(cluster_neurons.auditory{3},:)))
xlim([-200 2750]);  box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,2,7)
plot(-1000:5000, nanmean(seq_sdf_out(cluster_neurons.auditory{4},:)))
xlim([-200 2750]);  box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

% Frontal
subplot(4,2,2)
plot(-1000:5000, nanmean(seq_sdf_out(cluster_neurons.frontal{1},:)))
xlim([-200 2750]);  box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,2,4)
plot(-1000:5000, nanmean(seq_sdf_out(cluster_neurons.frontal{2},:)))
xlim([-200 2750]);  box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,2,6)
plot(-1000:5000, nanmean(seq_sdf_out(cluster_neurons.frontal{3},:)))
xlim([-200 2750]);  box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

subplot(4,2,8)
plot(-1000:5000, nanmean(seq_sdf_out(cluster_neurons.frontal{4},:)))
xlim([-200 2750]);  box off
vline(sound_onset_ms, 'k-')
vline([413	976	1539	2102	2665], 'k--')

%% Run cluster average autocorrelation

for g = 1:4
    % Compute population-averaged autocorrelation for auditory and frontal neurons
    [acorr_aud_pos(g,:), lags] = xcorr(nanmean(seq_sdf_out(cluster_neurons.auditory{g},acorr_time)), 'coeff');
    [acorr_frontal_pos(g,:), lags] = xcorr(nanmean(seq_sdf_out(cluster_neurons.frontal{g},acorr_time)), 'coeff');
end

figuren('Renderer','painters','Position',[100 100 250 500]);
subplot(2,1,1)
plot(lags, acorr_aud_pos)
xlim([-1000 1000]); ylim([-0.4 1]); box off; axis square

subplot(2,1,2)
plot(lags, acorr_frontal_pos)
xlim([-1000 1000]); ylim([-0.4 1]); box off; axis square

%%

nboot = 1000;
nSample = 100;

for g = 1:4
    clear peak_autocorr
    fprintf('Processing cluster group %i \n', g)
    for boot_i = 1:nboot

        [acorr_aud_pos_boot] = xcorr(nanmean(seq_sdf_out(randsample(cluster_neurons.auditory{g}, nSample, true),acorr_time)), 'coeff');
        [acorr_frontal_pos_boot] = xcorr(nanmean(seq_sdf_out(randsample(cluster_neurons.frontal{g}, nSample, true),acorr_time)), 'coeff');

        clear peak_idx_aud peak_idx_frontal 
        [~,peak_idx_aud] = findpeaks(acorr_aud_pos_boot,'MinPeakProminence',0.25);
        [~,peak_idx_frontal] = findpeaks(acorr_frontal_pos_boot,'MinPeakProminence',0.25);

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

    cluster_autocorr{g} = peak_autocorr;
end

%%

peak_ms_detected = []; peak_label = []; peak_area = [];
for g = 1:4
    p_peaks_detected(g,:) = mean(~isnan(cluster_autocorr{g}))*100;
    peak_ms_detected = [peak_ms_detected; cluster_autocorr{g}];
    peak_label = [peak_label; repmat(group_labels(g), nboot,1)];
    peak_area = [peak_area; repmat({'Auditory','Frontal'}, nboot,1)];
end

figuren('Renderer','painters','Position',[100 100 300 400]);
bar(p_peaks_detected)

figure;
plot_peaks = gramm('x', [peak_area(:,1); peak_area(:,2)], 'y', [peak_ms_detected(:,1); peak_ms_detected(:,2)], 'color', [peak_label; peak_label]);
plot_peaks.stat_summary('geom',{'point','errorbar'});
plot_peaks.axe_property('YLim',[500 1000]);
plot_peaks.geom_hline('yintercept',[1 2] * 563);
plot_peaks.draw

%%

disp('Auditory - facilitated')
comp_bootstrap_fixedvalue(cluster_autocorr{1}(:,1), 563)

disp('Frontal - facilitated')
comp_bootstrap_fixedvalue(cluster_autocorr{1}(:,2), 563)

disp('Auditory - suppressed')
comp_bootstrap_fixedvalue(cluster_autocorr{2}(:,1), 563)

disp('Frontal - suppressed')
comp_bootstrap_fixedvalue(cluster_autocorr{2}(:,2), 563)

