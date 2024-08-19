%% Figure: plot example neuron raster, example neuron sdf, and population sdf
% //// Fig 1a - Auditory cortex | Positive
example_neuron = 1498;
example_neuron_data = load(fullfile(dirs.mat_data, spike_log.session{example_neuron}));

% Align spikes and generate SDF
aligntime = example_neuron_data.event_table.stimulusOnset_ms;
for row = 1:length(aligntime)
    aligntime(row,:) = aligntime(row,:); + session_audio_latency{find(strcmp(spike_log.session{example_neuron}, ephysLog.session))}(row);
end

ops.sdf_filter = 'Gauss';
[sdf, raster] = get_spikes_aligned(example_neuron_data.spikes,aligntime,ops);
clear spk_figure_idv plot_ind_seq plot_pop_seq

% Example raster
spk_figure_idv(1,1)=gramm('x',raster.(spike_log.unitDSP{example_neuron}));
spk_figure_idv(1,1).geom_raster('geom',{'line'});
spk_figure_idv(1,1).axe_property('XLim',[-500 3000],'YLim',[0 size(raster.(spike_log.unitDSP{example_neuron}),1)]);
spk_figure_idv(1,1).geom_vline('xintercept',[0 563 1126 1689 2252]);

% Example SDF
for row = 1:size(sdf.(spike_log.unitDSP{example_neuron}),1)
    plot_ind_seq(row,:) = smooth(sdf.(spike_log.unitDSP{example_neuron})(row,:), 50);
end
spk_figure_idv(2,1)=gramm('x',ops.timewin,'y',plot_ind_seq);
spk_figure_idv(2,1).stat_summary();
spk_figure_idv(2,1).axe_property('XLim',[-500 3000],'YLim',[0 50]);

% Population SDF
for row = 1:length(sig_pos_units_auditory)
    plot_pop_seq(row,:) = smooth(norm_fr_sequence(sig_pos_units_auditory(row),:), 50);
end
spk_figure_idv(3,1)=gramm('x',ops.timewin,'y',plot_pop_seq);
spk_figure_idv(3,1).stat_summary();
spk_figure_idv(3,1).axe_property('XLim',[-500 3000],'YLim',[-0.5 1.5]);

% Figure layout
spk_figure_idv(1,1).set_layout_options...
    ('Position',[0.2 0.8 0.8 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);
spk_figure_idv(2,1).set_layout_options...
    ('Position',[0.2 0.45 0.8 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);
spk_figure_idv(3,1).set_layout_options...
    ('Position',[0.2 0.1 0.8 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);
% Draw figure
fig = figure('Renderer', 'painters', 'Position', [100 100 300 800]);
spk_figure_idv.draw();