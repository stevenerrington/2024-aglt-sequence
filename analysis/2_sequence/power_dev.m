% Parameters for pwelch
window = 500; % Length of each segment
noverlap = 250; % Number of overlapping samples
nfft = 5000; % Number of FFT points

clear power_*
parfor neuron_i = 1:size(spike_log,1)

    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

    sdf_in = load(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-laminar\data\spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
    event_table_in = load(fullfile(dirs.mat_data,[spike_log.session{neuron_i} '.mat']),'event_table');

    sdf_out = []; sdf_bl_out = [];
    for trial_i = 1:size(event_table_in.event_table,1)
        sdf_trial = []; sdf_trial = smooth(sdf_in.sdf.sequenceOnset(trial_i,1000+[0:3000]),100)';

        sdf_out(trial_i,:) = smooth(sdf_in.sdf.sequenceOnset(trial_i,:),100)';
        sdf_bl_out(trial_i,:) = smooth(sdf_in.sdf.reward(trial_i,:),100)';
    end

    neuron_seq_concat = []; neuron_bl_concat = [];

    for seq_i = 1:16
        neuron_seq_concat = [neuron_seq_concat, nanmean(sdf_out(event_table_in.event_table.cond_value == seq_i, 1000+[0:2750]))];
        neuron_bl_concat = [neuron_bl_concat, nanmean(sdf_bl_out(event_table_in.event_table.cond_value == seq_i, 1000+[500:3250]))];
    end

    neuron_seq_concat = neuron_seq_concat(~isnan(neuron_seq_concat));
    neuron_bl_concat = neuron_bl_concat(~isnan(neuron_bl_concat));

    [power_sequence(:,neuron_i), f(:,neuron_i)] = pwelch(neuron_seq_concat, window, noverlap, nfft, 1000, 'psd');
    [power_bl(:,neuron_i), ~] = pwelch(neuron_bl_concat, window, noverlap, nfft, 1000, 'psd');
end

power_sequence = power_sequence';
power_bl = power_bl';

cell_class = 'modulated';




power_idx = find(f(:,1) > 1.5 & f(:,1) < 2);

figuren('Renderer', 'painters', 'Position', [100 100 700 300]);
subplot(1,2,1);hold on
histogram(sum(power_sequence(neuron_class.auditory.(cell_class),power_idx),2)./sum(power_bl(neuron_class.auditory.(cell_class),power_idx),2),0:0.1:2,'LineStyle','None')
histogram(sum(power_sequence(neuron_class.frontal.(cell_class),power_idx),2)./sum(power_bl(neuron_class.frontal.(cell_class),power_idx),2),0:0.1:2,'LineStyle','None')
vline(1,'k--')
xlabel('Power ratio'); ylabel('Neurons')

subplot(1,2,2);hold on
scatter(sum(power_bl(neuron_class.auditory.(cell_class),power_idx),2),sum(power_sequence(neuron_class.auditory.(cell_class),power_idx),2),'filled')
scatter(sum(power_bl(neuron_class.frontal.(cell_class),power_idx),2),sum(power_sequence(neuron_class.frontal.(cell_class),power_idx),2),'filled')
plot(0:2000,0:2000,'k--');
xlabel('Baseline power'); ylabel('Sequence power')