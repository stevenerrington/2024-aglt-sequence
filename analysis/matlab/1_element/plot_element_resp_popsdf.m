%% Figure: plot example neuron raster, example neuron sdf, and population sdf
% //// Fig 1a - Auditory cortex | Positive
sound_list = {'A','C','D','F','G'};

xlim_vals = [-100 500];
ylim_vals = [0 60];
clear single_unit_fig

% Auditory
clear pop_
pop_neurons_in_aud = [aud_mod_neurons_pos; aud_mod_neurons_neg];
pop_label_aud = [repmat({'Positive'},length(aud_mod_neurons_pos),1); repmat({'Negative'},length(aud_mod_neurons_neg),1)];
pop_sdf_in_aud = [normFR_in.norm_fr_soundA(pop_neurons_in_aud,:); normFR_in.norm_fr_soundC(pop_neurons_in_aud,:); normFR_in.norm_fr_soundD(pop_neurons_in_aud,:); normFR_in.norm_fr_soundF(pop_neurons_in_aud,:); normFR_in.norm_fr_soundG(pop_neurons_in_aud,:)];

for i = 1:size(pop_sdf_in_aud,1)
    pop_sdf_in_aud(i,:) = smooth(pop_sdf_in_aud(i,:),1);
end

% Frontal
pop_neurons_in_frontal = [frontal_mod_neurons_pos; frontal_mod_neurons_neg];
pop_label_frontal = [repmat({'Positive'},length(frontal_mod_neurons_pos),1); repmat({'Negative'},length(frontal_mod_neurons_neg),1)];
pop_sdf_in_frontal = [normFR_in.norm_fr_soundA(pop_neurons_in_frontal,:); normFR_in.norm_fr_soundC(pop_neurons_in_frontal,:); normFR_in.norm_fr_soundD(pop_neurons_in_frontal,:); normFR_in.norm_fr_soundF(pop_neurons_in_frontal,:); normFR_in.norm_fr_soundG(pop_neurons_in_frontal,:)];

for i = 1:size(pop_sdf_in_frontal,1)
    pop_sdf_in_frontal(i,:) = smooth(pop_sdf_in_frontal(i,:),1);
end


pop_sound_cond_aud = [repmat({'A'}, length(pop_neurons_in_aud),1); repmat({'C'}, length(pop_neurons_in_aud),1); repmat({'D'}, length(pop_neurons_in_aud),1); repmat({'F'}, length(pop_neurons_in_aud),1); repmat({'G'}, length(pop_neurons_in_aud),1)];
pop_sound_cond_frontal = [repmat({'A'}, length(pop_neurons_in_frontal),1); repmat({'C'}, length(pop_neurons_in_frontal),1); repmat({'D'}, length(pop_neurons_in_frontal),1); repmat({'F'}, length(pop_neurons_in_frontal),1); repmat({'G'}, length(pop_neurons_in_frontal),1)];

single_unit_fig(1,1)=gramm('x',ops.sound_sdf_window,'y',pop_sdf_in_aud,'color',pop_sound_cond_aud,'row',repmat(pop_label_aud,5,1));
single_unit_fig(1,1).stat_summary();
single_unit_fig(1,1).axe_property('XLim',xlim_vals,'YLim',[-1 1.5]);


single_unit_fig(1,2)=gramm('x',ops.sound_sdf_window,'y',pop_sdf_in_frontal,'color',pop_sound_cond_frontal,'row',repmat(pop_label_frontal,5,1));
single_unit_fig(1,2).stat_summary();
single_unit_fig(1,2).axe_property('XLim',xlim_vals,'YLim',[-1 1.5]);


% Figure setup    ////////////////////////////////////////////////////////

single_unit_fig(1,1).geom_vline('xintercept',0,'style','k-');
single_unit_fig(1,1).set_names('x','Time from event (ms)','y','FR (spk/sec)');
single_unit_fig(1,2).geom_vline('xintercept',0,'style','k-');
single_unit_fig(1,2).set_names('x','Time from event (ms)','y','FR (spk/sec)');


figure('Renderer', 'painters', 'Position', [100 100 800 600]);
single_unit_fig.draw();