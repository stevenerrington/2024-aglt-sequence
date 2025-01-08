
%%
sound_times = [0, 563, 1126, 1689, 2252];

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
    nonviol_sdf = sdf_in.sdf.sequenceOnset(strcmp(event_table_in.event_table.cond_label, 'nonviol'), :);

    % Calculate baseline firing rate mean and standard deviation over a pre-stimulus period (-200 ms to 0 ms)
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, [1000+[-200:0]])));
    baseline_fr_std = nanstd(nanmean(nonviol_sdf(:, [1000+[-200:0]])));

    % Compute and smooth normalized SDF for the current neuron across time points
    seq_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 50);
end



%%

acorr_time = 1000+[0:2665];

[acorr_aud_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.auditory.glm_pos,acorr_time)), 'coeff');
[acorr_aud_neg, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.auditory.glm_neg,acorr_time)), 'coeff');
[acorr_frontal_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.frontal.glm_pos,acorr_time)), 'coeff');
[acorr_frontal_neg, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.frontal.glm_neg,acorr_time)), 'coeff');

[acorr_frontal_nonmod, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.frontal.nonmodulated,acorr_time)), 'coeff');
[acorr_auditory_nonmod, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.auditory.nonmodulated,acorr_time)), 'coeff');

figuren('Renderer', 'painters', 'Position', [100 100 500 700]); hold on
subplot(2,1,1); hold on; box off
plot(ops.timewin,nanmean(seq_sdf_out(neuron_class.auditory.glm_pos,:)))
plot(ops.timewin,nanmean(seq_sdf_out(neuron_class.auditory.glm_neg,:)))
plot(ops.timewin,nanmean(seq_sdf_out(neuron_class.frontal.glm_pos,:)))
plot(ops.timewin,nanmean(seq_sdf_out(neuron_class.frontal.glm_neg,:)))
xlim([-250 2750]); vline(sound_times,'k--')
legend({'Auditory: positive','Auditory: negative','Frontal: positive','Frontal: negative'})

subplot(2,1,2); hold on; box off
plot(lags,acorr_aud_pos)
plot(lags,acorr_aud_neg)
plot(lags,acorr_frontal_pos)
plot(lags,acorr_frontal_neg)
xlim([-1250 1250]); vline([-1126 -563 0 563 1126],'k--')
legend({'Auditory: positive','Auditory: negative','Frontal: positive','Frontal: negative'})


count =  0;

for area = {'frontal','auditory'}
    for type = {'glm_pos','glm_neg'}
        bootstrap_lags = []; mean_autocorr_r = [];

        for bootstrap_i = 1:1000
            % is it legal to resample from the same pool, or is it double dipping?

            clear acorr lags Xpk peak_lags
            [acorr, lags] = xcorr(nanmean(seq_sdf_out(randsample(neuron_class.(area{1}).(type{1}), 10, false),acorr_time)), 'coeff');
            lag_idx = find(lags > 500 & lags < 650);

            [~,Xpk,~,~] = findpeaks(acorr,'MinPeakProminence',0.2);
            peak_lags = lags(Xpk);

            try
                bootstrap_lags(bootstrap_i,1) = peak_lags(find(peak_lags > 0 , 1, 'first'));
            catch
                bootstrap_lags(bootstrap_i,1) = NaN;
            end

            mean_autocorr_r(bootstrap_i) = nanmean(acorr(lag_idx));

        end

        count = count + 1;
        label{count} = [int2str(count) '_' area{1} '_' type{1}];
        acorr_r{count} = mean_autocorr_r;
        acorr_peaks{count} = bootstrap_lags;


    end
end

figuren; 
subplot(2,1,1); hold on
histogram(acorr_r{1},-1:0.02:1,'LineStyle','none')
histogram(acorr_r{2},-1:0.02:1,'LineStyle','none')
histogram(acorr_r{3},-1:0.02:1,'LineStyle','none')
histogram(acorr_r{4},-1:0.02:1,'LineStyle','none')
legend(label, 'Interpreter', 'none')

subplot(2,1,2); hold on
histogram(acorr_peaks{1},0:10:800,'LineStyle','none')
histogram(acorr_peaks{2},0:10:800,'LineStyle','none')
histogram(acorr_peaks{3},0:10:800,'LineStyle','none')
histogram(acorr_peaks{4},0:10:800,'LineStyle','none')

mode(acorr_r{1}(~isnan(acorr_r{1})))
mode(acorr_r{2}(~isnan(acorr_r{2})))
mode(acorr_r{3}(~isnan(acorr_r{3})))
mode(acorr_r{4}(~isnan(acorr_r{4})))


mean(~isnan(acorr_peaks{1}))*100
mean(~isnan(acorr_peaks{2}))*100
mean(~isnan(acorr_peaks{3}))*100
mean(~isnan(acorr_peaks{4}))*100


