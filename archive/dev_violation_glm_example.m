%%

neuron_i = 1535;

% Load the spike data for the current neuron
sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

% Load the event table for the current session
event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

cond_inc = [];
cond_inc = [1 5 14 2 6 15 3 7 13 4 8 16];


for trial_i = 1:size(sdf_in.sdf.violation,1)
    sdf_in.sdf.violation(trial_i,:) = smooth(sdf_in.sdf.violation(trial_i,:),50);
end


clear viol_sdf_zscore nonviol_sdf_zscore
for neuron_i = 1:2298
    viol_sdf_zscore(neuron_i,:) = smooth(viol_sdf(neuron_i,:),50)';
    nonviol_sdf_zscore(neuron_i,:) = smooth(nonviol_sdf(neuron_i,:),50)';
end


clear spk_figure_idv

% Example raster
spk_figure_idv(1,1)=gramm('x',sdf_in.raster.violation,'color',event_table_in.event_table.cond_label,'subset',ismember(event_table_in.event_table.cond_value,cond_inc));
spk_figure_idv(1,1).geom_raster('geom',{'line'});
spk_figure_idv(1,1).axe_property('XLim',[-250 800]);
spk_figure_idv(1,1).geom_vline('xintercept',[0 413 563]);

% Example SDF
spk_figure_idv(2,1)=gramm('x',ops.timewin,'y',sdf_in.sdf.violation,'color',event_table_in.event_table.cond_label,'subset',ismember(event_table_in.event_table.cond_value,cond_inc));
spk_figure_idv(2,1).stat_summary();
spk_figure_idv(2,1).axe_property('XLim',[-250 800],'YLim',[0 15]);
spk_figure_idv(2,1).geom_vline('xintercept',[0 413 563]);

% Population
spk_figure_idv(3,1)=gramm('x',ops.sound_sdf_window,'y',[viol_sdf_zscore; nonviol_sdf_zscore], 'color', [repmat({'viol'},size(viol_sdf_zscore,1),1); repmat({'nonviol'},size(nonviol_sdf_zscore,1),1)], 'subset', repmat(ismember(1:2298,auditory_viol_neurons_neg)',2,1));
spk_figure_idv(3,1).stat_summary();
spk_figure_idv(3,1).axe_property('XLim',[-250 800],'YLim',[-2 1]);
spk_figure_idv(3,1).geom_vline('xintercept',[0 413 563]);


% Figure layout
spk_figure_idv(1,1).set_layout_options...
    ('Position',[0.1 0.8 0.8 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

spk_figure_idv(2,1).set_layout_options...
    ('Position',[0.1 0.45 0.8 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

spk_figure_idv(3,1).set_layout_options...
    ('Position',[0.1 0.1 0.8 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);


% Draw figure
fig = figure('Renderer', 'painters', 'Position', [100 100 400 800]);
spk_figure_idv.draw();
