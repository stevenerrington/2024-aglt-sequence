clear glm_seq* encoding_seq_*

% Parallel loop to perform GLM analysis for each neuron
parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

    % Perform GLM analysis on sound-aligned SDF data for the current neuron
    [glm_seq_output{neuron_i},...
        glm_seq_encoding_flag_ord(neuron_i,:), glm_seq_encoding_flag_rel(neuron_i,:),...
        encoding_seq_beta_ord(neuron_i,:), encoding_seq_beta_rel(neuron_i,:)] =...
        glm_position_relative(sdf_soundAlign_data{neuron_i});
end

for neuron_i = 1:size(spike_log,1)
    abs_pos_rvalue(neuron_i,1) = nanmax(glm_seq_output{neuron_i}.ordinal.var_exp(16:56));
    rel_pos_rvalue(neuron_i,1) = nanmax(glm_seq_output{neuron_i}.relative.var_exp(16:56));
end

%% /////////////////////////////////////////////////////////////////

figuren('Renderer', 'painters', 'Position', [100 100 300 400]); 
subplot(2,1,1); hold on; box off
histogram(abs_pos_rvalue(neuron_class.auditory.modulated),-0.1:0.01:0.3,'LineStyle','None')
histogram(rel_pos_rvalue(neuron_class.auditory.modulated),-0.1:0.01:0.3,'LineStyle','None')

subplot(2,1,2); hold on; box off
histogram(abs_pos_rvalue(neuron_class.frontal.modulated),-0.1:0.01:0.3,'LineStyle','None')
histogram(rel_pos_rvalue(neuron_class.frontal.modulated),-0.1:0.01:0.3,'LineStyle','None')

ttest(abs_pos_rvalue(neuron_class.auditory.modulated), rel_pos_rvalue(neuron_class.auditory.modulated))
ttest(abs_pos_rvalue(neuron_class.frontal.modulated), rel_pos_rvalue(neuron_class.frontal.modulated))

%% /////////////////////////////////////////////////////////////////

sig_relative_enc_neurons = sum(glm_seq_encoding_flag_rel,2) > 0;
sig_absolute_enc_neurons = sum(glm_seq_encoding_flag_ord,2) > 0;

sig_lookup_table = [sig_relative_enc_neurons, sig_absolute_enc_neurons];

disp(['Auditory | Both: ' int2str(sum(sum(sig_lookup_table(neuron_class.auditory.modulated,:),2) == 2)) ])
disp(['Auditory | Rel only: ' int2str(sum(sig_lookup_table(neuron_class.auditory.modulated,1) == 1 & sig_lookup_table(neuron_class.auditory.modulated,2) == 0 )) ])
disp(['Auditory | Abs only: ' int2str(sum(sig_lookup_table(neuron_class.auditory.modulated,1) == 0 & sig_lookup_table(neuron_class.auditory.modulated,2) == 1 )) ])
disp(['Auditory | Neither: ' int2str(sum(sig_lookup_table(neuron_class.auditory.modulated,1) == 0 & sig_lookup_table(neuron_class.auditory.modulated,2) == 0 )) ])


disp(['Frontal | Both: ' int2str(sum(sum(sig_lookup_table(neuron_class.frontal.modulated,:),2) == 2)) ])
disp(['Frontal | Rel only: ' int2str(sum(sig_lookup_table(neuron_class.frontal.modulated,1) == 1 & sig_lookup_table(neuron_class.frontal.modulated,2) == 0 )) ])
disp(['Frontal | Abs only: ' int2str(sum(sig_lookup_table(neuron_class.frontal.modulated,1) == 0 & sig_lookup_table(neuron_class.frontal.modulated,2) == 1 )) ])
disp(['Frontal | Neither: ' int2str(sum(sig_lookup_table(neuron_class.frontal.modulated,1) == 0 & sig_lookup_table(neuron_class.frontal.modulated,2) == 0 )) ])


%% ////////////////////////////////////////////////////////////
neuron_i = 1498;

trials_in = []; trials_in = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C');
sdf_in = []; sdf_in = sdf_soundAlign_data{neuron_i}(trials_in,1);
sdf_out = []; average_fr = [];

for trial_i = 1:size(sdf_in,1)
    sdf_out(trial_i,:) = smooth(sdf_in{trial_i},50)';
    average_fr(trial_i,:) = nanmean(sdf_out(trial_i,200+[0:200]));
end

clear pos_rel_sdf_figure
pos_rel_sdf_figure(1, 1) = gramm('x', ops.sound_sdf_window, 'y', sdf_out, 'color', sdf_soundAlign_data{neuron_i}(trials_in,4));
pos_rel_sdf_figure(1, 1).stat_summary();
pos_rel_sdf_figure(1, 1).axe_property('XLim', [-100 500], 'YLim', [0 30]);
pos_rel_sdf_figure(1, 1).geom_vline('XIntercept',0);

pos_rel_sdf_figure(1, 2) = gramm('x', sdf_soundAlign_data{neuron_i}(trials_in,4), 'y', average_fr, 'color', sdf_soundAlign_data{neuron_i}(trials_in,4));
pos_rel_sdf_figure(1, 2).stat_summary('geom',{'lines'});
pos_rel_sdf_figure(1, 2).axe_property('YLim', [5 15]);

pos_rel_sdf_figure(1, 3) = gramm('x', ops.sound_sdf_window, 'y', sdf_out, 'color', sdf_soundAlign_data{neuron_i}(trials_in,6));
pos_rel_sdf_figure(1, 3).stat_summary();
pos_rel_sdf_figure(1, 3).axe_property('XLim', [-100 500], 'YLim', [0 30]);
pos_rel_sdf_figure(1, 3).geom_vline('XIntercept',0);

pos_rel_sdf_figure(1, 4) = gramm('x', sdf_soundAlign_data{neuron_i}(trials_in,6), 'y', average_fr, 'color', sdf_soundAlign_data{neuron_i}(trials_in,6));
pos_rel_sdf_figure(1, 4).stat_summary('geom',{'lines'});
pos_rel_sdf_figure(1, 4).axe_property('YLim', [5 15]);

pos_rel_sdf_figure(1,1).set_layout_options...
    ('Position',[0.1 0.15 0.35 0.7],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

pos_rel_sdf_figure(1,2).set_layout_options...
    ('Position',[0.35 0.7 0.1 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

pos_rel_sdf_figure(1,3).set_layout_options...
    ('Position',[0.55 0.15 0.35 0.7],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

pos_rel_sdf_figure(1,4).set_layout_options...
    ('Position',[0.8 0.7 0.1 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);


figure('Renderer', 'painters', 'Position', [100 100 1000 400]);
pos_rel_sdf_figure.draw'

