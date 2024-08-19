%% Figure: plot example neuron raster, example neuron sdf, and population sdf
% //// Fig 1a - Auditory cortex | Positive
sound_list = {'A','C','D','F','G'};

xlim_vals = [-100 500];
ylim_vals = [0 60];
clear single_unit_fig

% Auditory modulation //////////////////////////////////////////////////////
clear example_*
example_valid_idx = strcmp(string(sdf_soundAlign_data{aud_neuron_example}(:,5)) ,'nonviol') &  ~strcmp(string(sdf_soundAlign_data{aud_neuron_example}(:,3)) ,'Baseline');
example_sdf_in = cell2mat(sdf_soundAlign_data{aud_neuron_example}(example_valid_idx,1));

for i = 1:size(example_sdf_in,1)
    example_sdf_in(i,:) = smooth(example_sdf_in(i,:),50);
end

example_raster_in = sdf_soundAlign_data{aud_neuron_example}(example_valid_idx,2);
example_sound_cond = sdf_soundAlign_data{aud_neuron_example}(example_valid_idx,3);
example_order_cond = sdf_soundAlign_data{aud_neuron_example}(example_valid_idx,4);

clear pop_
pop_neurons_in = aud_mod_neurons;
pop_sdf_in = [norm_fr_soundA(pop_neurons_in,:); norm_fr_soundC(pop_neurons_in,:); norm_fr_soundD(pop_neurons_in,:); norm_fr_soundF(pop_neurons_in,:); norm_fr_soundG(pop_neurons_in,:)];

for i = 1:size(pop_sdf_in,1)
    pop_sdf_in(i,:) = smooth(pop_sdf_in(i,:),50);
end

pop_sound_cond = [repmat({'A'}, length(pop_neurons_in),1); repmat({'C'}, length(pop_neurons_in),1); repmat({'D'}, length(pop_neurons_in),1); repmat({'F'}, length(pop_neurons_in),1); repmat({'G'}, length(pop_neurons_in),1)];

single_unit_fig(1,1)=gramm('x',example_raster_in,'color',example_sound_cond);
single_unit_fig(1,1).geom_raster('geom',{'point'});
single_unit_fig(1,1).axe_property('XLim',xlim_vals);

single_unit_fig(2,1)=gramm('x',ops.sound_sdf_window,'y',example_sdf_in,'color',example_sound_cond);
single_unit_fig(2,1).stat_summary();
single_unit_fig(2,1).axe_property('XLim',xlim_vals,'YLim',[0 50]);

single_unit_fig(3,1)=gramm('x',ops.sound_sdf_window,'y',pop_sdf_in,'color',pop_sound_cond);
single_unit_fig(3,1).stat_summary();
single_unit_fig(3,1).axe_property('XLim',xlim_vals,'YLim',[-2 3]);

% Frontal modulation //////////////////////////////////////////////////////
clear example_*
example_valid_idx = strcmp(string(sdf_soundAlign_data{frontal_neuron_example}(:,5)) ,'nonviol') &  ~strcmp(string(sdf_soundAlign_data{frontal_neuron_example}(:,3)) ,'Baseline');
example_sdf_in = cell2mat(sdf_soundAlign_data{frontal_neuron_example}(example_valid_idx,1));

for i = 1:size(example_sdf_in,1)
    example_sdf_in(i,:) = smooth(example_sdf_in(i,:),50);
end

example_raster_in = sdf_soundAlign_data{frontal_neuron_example}(example_valid_idx,2);
example_sound_cond = sdf_soundAlign_data{frontal_neuron_example}(example_valid_idx,3);
example_order_cond = sdf_soundAlign_data{frontal_neuron_example}(example_valid_idx,4);

clear pop_
pop_neurons_in = frontal_mod_neurons;
pop_sdf_in = [norm_fr_soundA(pop_neurons_in,:); norm_fr_soundC(pop_neurons_in,:); norm_fr_soundD(pop_neurons_in,:); norm_fr_soundF(pop_neurons_in,:); norm_fr_soundG(pop_neurons_in,:)];

for i = 1:size(pop_sdf_in,1)
    pop_sdf_in(i,:) = smooth(pop_sdf_in(i,:),50);
end

pop_sound_cond = [repmat({'A'}, length(pop_neurons_in),1); repmat({'C'}, length(pop_neurons_in),1); repmat({'D'}, length(pop_neurons_in),1); repmat({'F'}, length(pop_neurons_in),1); repmat({'G'}, length(pop_neurons_in),1)];

single_unit_fig(1,2)=gramm('x',example_raster_in,'color',example_sound_cond);
single_unit_fig(1,2).geom_raster('geom',{'point'});
single_unit_fig(1,2).axe_property('XLim',xlim_vals);

single_unit_fig(2,2)=gramm('x',ops.sound_sdf_window,'y',example_sdf_in,'color',example_sound_cond);
single_unit_fig(2,2).stat_summary();
single_unit_fig(2,2).axe_property('XLim',xlim_vals,'YLim',[20 70]);

single_unit_fig(3,2)=gramm('x',ops.sound_sdf_window,'y',pop_sdf_in,'color',pop_sound_cond);
single_unit_fig(3,2).stat_summary();
single_unit_fig(3,2).axe_property('XLim',xlim_vals,'YLim',[-2 3]);


% Figure setup    ////////////////////////////////////////////////////////
single_unit_fig(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
single_unit_fig(2,1).axe_property('XTick',[],'XColor',[1 1 1]);
single_unit_fig(1,2).axe_property('XTick',[],'XColor',[1 1 1]);
single_unit_fig(2,2).axe_property('XTick',[],'XColor',[1 1 1]);


single_unit_fig(1,1).geom_vline('xintercept',0,'style','k-');
single_unit_fig(2,1).geom_vline('xintercept',0,'style','k-');
single_unit_fig(3,1).geom_vline('xintercept',0,'style','k-');

single_unit_fig(1,2).geom_vline('xintercept',0,'style','k-');
single_unit_fig(2,2).geom_vline('xintercept',0,'style','k-');
single_unit_fig(3,2).geom_vline('xintercept',0,'style','k-');


single_unit_fig(1,1).set_names('y','Trials');
single_unit_fig(2,1).set_names('x','','y','FR (spk/sec)');
single_unit_fig(3,1).set_names('x','Time from event (ms)','y','FR (spk/sec)');

single_unit_fig(1,2).set_names('y','Trials');
single_unit_fig(2,2).set_names('x','','y','FR (spk/sec)');
single_unit_fig(3,2).set_names('x','Time from event (ms)','y','FR (spk/sec)');


% Figure layout ////////////////////////////////////////////////////////////
single_unit_fig(1,1).set_layout_options...
    ('Position',[0.15 0.8 0.35 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(2,1).set_layout_options...
    ('Position',[0.15 0.45 0.35 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(3,1).set_layout_options...
    ('Position',[0.15 0.1 0.35 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(1,2).set_layout_options...
    ('Position',[0.6 0.8 0.35 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(2,2).set_layout_options...
    ('Position',[0.6 0.45 0.35 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(3,2).set_layout_options...
    ('Position',[0.6 0.1 0.35 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

figure('Renderer', 'painters', 'Position', [100 100 800 600]);
single_unit_fig.draw();