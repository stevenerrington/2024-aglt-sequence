
example_neuron_aud = 1286;
example_neuron_frontal = 1182;

example_neuron_in = example_neuron_frontal;

example_neuron_data = load(fullfile(dirs.mat_data, spike_log.session{example_neuron_in}));
ops.plot_ch = spike_log.unitDSP(example_neuron_in);

aligntime = example_neuron_data.event_table.stimulusOnset_ms;
ops.sdf_filter = 'Gauss';
[sdf, raster] = get_spikes_aligned(example_neuron_data.spikes,aligntime,ops);


% Example raster
spk_figure_idv(1,1)=gramm('x',raster.(spike_log.unitDSP{example_neuron_in}));
spk_figure_idv(1,1).geom_raster('geom',{'line'});
spk_figure_idv(1,1).axe_property('XLim',[-500 3500],'YLim',[0 size(raster.(spike_log.unitDSP{example_neuron_in}),1)]);
spk_figure_idv(1,1).geom_vline('xintercept',[0 563 1126 1689 2252]);

% Example SDF
spk_figure_idv(2,1)=gramm('x',ops.timewin,'y',sdf.(spike_log.unitDSP{example_neuron_in}));
spk_figure_idv(2,1).stat_summary();
spk_figure_idv(2,1).axe_property('XLim',[-500 3500],'YLim',[0 12]);
spk_figure_idv(2,1).geom_vline('xintercept',[0 563 1126 1689 2252]);


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



% Draw figure
fig = figure('Renderer', 'painters', 'Position', [100 100 600 800]);
spk_figure_idv.draw();




%%
example_neuron_aud = 1286;
example_neuron_frontal = 1182;
example_neuron_in = example_neuron_frontal;

example_neuron_data = load(fullfile(dirs.mat_data, spike_log.session{example_neuron_in}));
event_table = example_neuron_data.event_table;
% Get violation alignment
for trial_i = 1:size(event_table,1)
    switch event_table.cond_value(trial_i)
        case {3 7 14}
            event_table.violation_ms(trial_i) =  event_table.stimulusOnset_ms(trial_i) + 1127;
        case {4 8 15}
            event_table.violation_ms(trial_i) =  event_table.stimulusOnset_ms(trial_i) + 2253;
        case {1 5 13}
            event_table.violation_ms(trial_i) =  event_table.stimulusOnset_ms(trial_i) + 1127;
        case {2 6 16}
            event_table.violation_ms(trial_i) =  event_table.stimulusOnset_ms(trial_i) + 2253;
        case {9 10 11 12}
            event_table.violation_ms(trial_i) =  NaN;
    end
end

ops.plot_ch = spike_log.unitDSP(example_neuron_in);
aligntime = event_table.violation_ms;
ops.sdf_filter = 'Gauss';
ops.timewin = [-1000:5000];
[sdf, raster] = get_spikes_aligned(example_neuron_data.spikes,aligntime,ops);

neuron_list = fieldnames(sdf);

ch = neuron_list{4}(end-2:end);
clear test

test(1,1)=gramm('x',raster.(['DSP' ch]),'color',event_table.cond_label);
test(1,1).geom_raster('geom',{'line'});
test(1,1).axe_property('XLim',[-200 600]);

test(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]),'color',event_table.cond_label);
test(2,1).stat_summary();
test(2,1).axe_property('XLim',[-200 600],'YLim',[0 20]);

test(1,1).axe_property('XTick',[],'XColor',[1 1 1]);

test(1,1).geom_vline('xintercept',0,'style','k-');
test(2,1).geom_vline('xintercept',0,'style','k-');

test(1,1).set_names('y','Trials');
test(2,1).set_names('x','Time from event (ms)','y','FR (spk/sec)');

test(1,1).set_layout_options...
    ('Position',[0.1 0.8 0.8 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

test(2,1).set_layout_options...
    ('Position',[0.1 0.1 0.8 0.6],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);


test.set_title(['DSP' ch]);
figure('Renderer', 'painters', 'Position', [100 100 700 600]);
test.draw();
