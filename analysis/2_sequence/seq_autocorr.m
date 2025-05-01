rng(1,"twister");

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
    seq_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);
end


%% Run autocorrelation

acorr_time = 1000+[0:2665];

[acorr_aud_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.auditory.all,acorr_time)), 'coeff');
[acorr_frontal_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.frontal.all,acorr_time)), 'coeff');


[acorr_clu1_auditory, lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.clu1,neuron_class.auditory.all),acorr_time)), 'coeff');
[acorr_clu2_auditory, lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.clu2,neuron_class.auditory.all),acorr_time)), 'coeff');
[acorr_clu3_auditory, lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.clu3,neuron_class.auditory.all),acorr_time)), 'coeff');
[acorr_clu4_auditory, lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.clu4,neuron_class.auditory.all),acorr_time)), 'coeff');

[acorr_clu1_frontal, lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.clu1,neuron_class.frontal.all),acorr_time)), 'coeff');
[acorr_clu2_frontal, lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.clu2,neuron_class.frontal.all),acorr_time)), 'coeff');
[acorr_clu3_frontal, lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.clu3,neuron_class.frontal.all),acorr_time)), 'coeff');
[acorr_clu4_frontal, lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.clu4,neuron_class.frontal.all),acorr_time)), 'coeff');

%% Figures
figure_xlim = [-250 2665];
figuren('Renderer', 'painters', 'Position', [100 100 700 400]); hold on

plot_ylim_clu_i = {[-0.4 0.4], [-0.25 0.25];...
    [-0.2 0.8], [-0.2 0.8];...
    [-0.6 0.6], [-0.6 0.2];...
    [-0.4 0.4], [-0.4 0.1]};


for clu_i = 1:4
    subplot(2,2,clu_i); hold on
    yyaxis left
    plot(-1000:5000, nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.auditory.all),:)),'LineWidth',figure_linewidth,'Color',color_pal.(['clu' int2str(clu_i)]))
    ylim(plot_ylim_clu_i{clu_i,1})
    
    yyaxis right
    plot(-1000:5000, nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.frontal.all),:)),'LineWidth',figure_linewidth+1,'Color',color_pal.(['clu' int2str(clu_i)]))
    ylim(plot_ylim_clu_i{clu_i,2})
    
    
    vline([sound_onset_ms], 'k-')
    vline([sound_onset_ms]+413, 'k--')
    xlim(figure_xlim)
end

figuren('Renderer', 'painters', 'Position', [100 100 700 400]); hold on
subplot(1,2,1); hold on; box off
plot(lags,acorr_clu1_auditory,'color',color_pal.clu1,'LineWidth',1.5); 
plot(lags,acorr_clu2_auditory,'color',color_pal.clu2,'LineWidth',1.5); 
plot(lags,acorr_clu3_auditory,'color',color_pal.clu3,'LineWidth',1.5); 
plot(lags,acorr_clu4_auditory,'color',color_pal.clu4,'LineWidth',1.5); 
xlim([-1250 1250]); vline([-1126 -563 0 563 1126],'k--')
title('Auditory')
legend({'Facilitated','Ramping','Suppressed','Peri-onset'},'Location','southoutside','NumColumns', 2)

subplot(1,2,2); hold on; box off
plot(lags,acorr_clu1_frontal,'color',color_pal.clu1,'LineWidth',1.5); 
plot(lags,acorr_clu2_frontal,'color',color_pal.clu2,'LineWidth',1.5); 
plot(lags,acorr_clu3_frontal,'color',color_pal.clu3,'LineWidth',1.5); 
plot(lags,acorr_clu4_frontal,'color',color_pal.clu4,'LineWidth',1.5); 
xlim([-1250 1250]); vline([-1126 -563 0 563 1126],'k--')
title('Frontal')
legend({'Facilitated','Ramping','Suppressed','Peri-onset'},'Location','southoutside','NumColumns', 2)



%%
clear label acorr_r acorr_peaks
count =  0;

for area = 1:2
    for type = {'A'}
        bootstrap_lags = []; mean_autocorr_r = [];

        for bootstrap_i = 1:1000
            % is it legal to resample from the same pool, or is it double dipping?
            sample_pop_idx = []; 
            
            switch area
                case 1
                    area_label = 'aud';
                    sample_pop_idx = neuron_class.auditory.all;
                case 2
                    area_label = 'frontal';
                    sample_pop_idx = neuron_class.frontal.all;
           end


            clear acorr lags Xpk peak_lags
            [acorr, lags] = xcorr(nanmean(seq_sdf_out(randsample(sample_pop_idx, 10, false),acorr_time)), 'coeff');
            lag_idx = find(lags > 500 & lags < 650);

            [~,Xpk,~,~] = findpeaks(acorr,'MinPeakProminence',0.33);
            peak_lags = lags(Xpk);

            try
                bootstrap_lags(bootstrap_i,1) = peak_lags(find(peak_lags > 0 , 1, 'first'));
            catch
                bootstrap_lags(bootstrap_i,1) = NaN;
            end

            mean_autocorr_r(bootstrap_i) = nanmean(acorr(lag_idx));

        end

        count = count + 1;
        label{count} = [int2str(count) '_' area_label '_' type{1}];
        acorr_r{count} = mean_autocorr_r;
        acorr_peaks{count} = bootstrap_lags;


    end
end

[nanmean(acorr_peaks{1}) nanmean(acorr_peaks{2})]
[nanstd(acorr_peaks{1}) nanstd(acorr_peaks{2})]


[mean(~isnan(acorr_peaks{1}))*100, mean(~isnan(acorr_peaks{2}))*100]

[h, p, ci, stats] = vartest2(acorr_peaks{1}, acorr_peaks{2});


[h,p,ci,stats] = ttest(acorr_peaks{1}, 563);
[h,p,ci,stats] = ttest(acorr_peaks{2}, 563);


%%
count =  0;

for area = 1:2
    for type = {'clu1','clu2','clu3','clu4'}
        bootstrap_lags = []; mean_autocorr_r = [];

        for bootstrap_i = 1:1000
            % is it legal to resample from the same pool, or is it double dipping?
            sample_pop_idx = []; 
            
            switch area
                case 1
                    area_label = 'aud';
                    sample_pop_idx = intersect(auditory_neuron_idx, neuron_class.cluster_idx.(type{1}));
                case 2
                    area_label = 'frontal';
                    sample_pop_idx = intersect(frontal_neuron_idx, neuron_class.cluster_idx.(type{1}));
           end


            clear acorr lags Xpk peak_lags
            [acorr, lags] = xcorr(nanmean(seq_sdf_out(randsample(sample_pop_idx, 10, false),acorr_time)), 'coeff');
            lag_idx = find(lags > 500 & lags < 650);

            [~,Xpk,~,~] = findpeaks(acorr,'MinPeakProminence',0.33);
            peak_lags = lags(Xpk);

            try
                bootstrap_lags(bootstrap_i,1) = peak_lags(find(peak_lags > 0 , 1, 'first'));
            catch
                bootstrap_lags(bootstrap_i,1) = NaN;
            end

            mean_autocorr_r(bootstrap_i) = nanmean(acorr(lag_idx));

        end

        count = count + 1;
        label{count} = [int2str(count) '_' area_label '_' type{1}];
        acorr_r{count} = mean_autocorr_r;
        acorr_peaks{count} = bootstrap_lags;


    end
end


acorr_peak_plotData = [acorr_peaks{1}; acorr_peaks{2}; acorr_peaks{3}; acorr_peaks{4};...
    acorr_peaks{5}; acorr_peaks{6}; acorr_peaks{7}; acorr_peaks{8}];

acorr_peak_plotLabels = [repmat(label(1),1000,1); repmat(label(2),1000,1); repmat(label(3),1000,1); repmat(label(4),1000,1);...
    repmat(label(5),1000,1); repmat(label(6),1000,1); repmat(label(7),1000,1); repmat(label(8),1000,1)];

% Plot 4: Population beta weight plot
xcorr_peak_time(1, 1) = gramm('x', acorr_peak_plotLabels, 'y', acorr_peak_plotData, 'color', acorr_peak_plotLabels);
xcorr_peak_time(1, 1).geom_hline('yintercept',[563 1126])
xcorr_peak_time(1, 1).stat_boxplot('width',6)
xcorr_peak_time(1, 1).geom_jitter()
xcorr_peak_time(1, 1).axe_property('YLim',[0 1500])

figure('Renderer', 'painters', 'Position', [100 100 500 300]);
xcorr_peak_time.draw


[mean(~isnan(acorr_peaks{1}))*100, mean(~isnan(acorr_peaks{2}))*100, mean(~isnan(acorr_peaks{3}))*100, mean(~isnan(acorr_peaks{4}))*100]
[mean(~isnan(acorr_peaks{5}))*100, mean(~isnan(acorr_peaks{6}))*100, mean(~isnan(acorr_peaks{7}))*100, mean(~isnan(acorr_peaks{8}))*100]


[nanmean(acorr_peaks{8}) nanstd(acorr_peaks{8})]

