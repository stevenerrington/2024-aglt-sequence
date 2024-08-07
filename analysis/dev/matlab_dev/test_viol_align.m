
example_neuron_aud = 1286;
example_neuron_frontal = 1182;

example_neuron_in = example_neuron_aud;

example_neuron_data = load(fullfile(dirs.mat_data, spike_log.session{example_neuron_in}));
ops.plot_ch = spike_log.unitDSP(example_neuron_in);

event_table = example_neuron_data.event_table;

% Get comparison group label
for trl_i = 1:size(event_table,1)
    if ismember(event_table.cond_value(trl_i), [3 7 14])
        event_table.comparison_label{trl_i} = 'single_viol_cond_a';
    elseif ismember(event_table.cond_value(trl_i), [4 8 15])
        event_table.comparison_label{trl_i} = 'single_viol_cond_b';
    elseif ismember(event_table.cond_value(trl_i), [1 5 13])
        event_table.comparison_label{trl_i} = 'single_viol_cond_c';
    elseif ismember(event_table.cond_value(trl_i), [2 6 16])
        event_table.comparison_label{trl_i} = 'single_viol_cond_d';
    else
        event_table.comparison_label{trl_i} = 'na';
    end
end

for trl_i = 1:size(event_table,1)
    event_table.stimulusOnset_ms(trl_i) = event_table.stimulusOnset_ms(trl_i) + 300;
end

aligntime = event_table.stimulusOnset_ms;

ops.timewin = -1000:5000;
ops.sdf_filter = 'Gauss';

[sdf, raster] = get_spikes_aligned(example_neuron_data.spikes,aligntime,ops);
names = fieldnames( example_neuron_data.spikes.time );

violation_win = [-500:1000];
xlim_vals = [violation_win(1) violation_win(end)];
ylim_vals = [0 25];


comp_group_list = {'single_viol_cond_a', 'single_viol_cond_b', 'single_viol_cond_c', 'single_viol_cond_d'};
comp_offset_list = [1127 2253 1127 2253];


for ch_i = 1:30
    ch = names{ch_i}(end-2:end);

    bl_mean = nanmean(nanmean(sdf.(['DSP' ch])(:,[-500:0]+1000)));
    bl_sd = nanstd(nanmean(sdf.(['DSP' ch])(:,[-500:0]+1000)));
    sdf_in = (sdf.(['DSP' ch]) - bl_mean)./bl_sd;

    viol_sdf = [];
    nonviol_sdf = [];

    for comparison_idx = 1:length(comp_group_list)

        comparson_group = comp_group_list{comparison_idx};

        viol_trl_in = []; viol_trl_in = find(strcmp(event_table.comparison_label,comparson_group) & strcmp(event_table.cond_label,'viol'));
        nonviol_trl_in = []; nonviol_trl_in = find(strcmp(event_table.comparison_label,comparson_group) & strcmp(event_table.cond_label,'nonviol'));

        viol_sdf = [viol_sdf; sdf_in(viol_trl_in,abs(ops.timewin(1))+violation_win+comp_offset_list(comparison_idx))];
        nonviol_sdf = [nonviol_sdf; sdf_in(nonviol_trl_in,abs(ops.timewin(1))+violation_win+comp_offset_list(comparison_idx))];
    end



    % figuren('Renderer', 'painters', 'Position', [100 100 700 600]); hold on
    % plot(violation_win,nanmean(viol_sdf),'color',[1 0 0],'LineWidth',1.5)
    % plot(violation_win,nanmean(nonviol_sdf),'color',[0 1 0],'LineWidth',1.5)

    pop_av_viol(ch_i,:) = nanmean(viol_sdf);
    pop_av_nonviol(ch_i,:) = nanmean(nonviol_sdf);

end

for ch_i = 1:30
    figuren('Renderer', 'painters', 'Position', [100 100 700 600]); hold on
    plot(violation_win,smooth(pop_av_viol(ch_i,:),50),'color',[1 0 0],'LineWidth',1.5)
    plot(violation_win,smooth(pop_av_nonviol(ch_i,:),50),'color',[0 1 0],'LineWidth',1.5)
end