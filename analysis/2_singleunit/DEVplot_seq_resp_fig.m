% Load data
neuron_count = 0; lfp_count = 0; % Initialize counters for neurons and LFPs

% Loop through each session in the ephysLog
for session_i = 1:size(ephysLog,1)

    % Load the session data file
    datafile = ephysLog.session{session_i};
    load(fullfile(dirs.mat_data,datafile))
    fprintf('Session %i of %i | %s \n', session_i, size(ephysLog,1), datafile)

    % Find all neurons recorded in the current session
    neuron_list = spike_log.unitDSP(strcmp(spike_log.session,datafile));

    % Process each neuron in the session
    for neuron_i = 1:length(neuron_list)
        neuron_label = neuron_list{neuron_i}; % Get the neuron label
        neuron_count = neuron_count + 1; % Increment the neuron count
        sound_sdf = {}; % Initialize sound SDF cell array
        count = 0; % Initialize counter for sound SDF

        load(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data\spike',[datafile '_' neuron_label]))


        baseline_fr_mu_trialStart = nanmean(nanmean(sdf.trialStart(:,1000+[-200:0])));
        baseline_fr_std_trialStart = nanstd(nanmean(sdf.trialStart(:,1000+[-200:0])));
        
        baseline_fr_mu_seq= nanmean(nanmean(sdf.sequenceOnset(:,1000+[-200:0])));
        baseline_fr_std_seq = nanstd(nanmean(sdf.sequenceOnset(:,1000+[-200:0])));
        
        baseline_fr_mu_reward = nanmean(nanmean(sdf.reward(:,1000+[-200:0])));
        baseline_fr_std_reward = nanstd(nanmean(sdf.reward(:,1000+[-200:0])));


        sdf_nonviol_trialStart(neuron_count,:) = (nanmean(sdf.trialStart(strcmp(event_table.cond_label,'nonviol'),:))-baseline_fr_mu_trialStart)./baseline_fr_std_trialStart;
        sdf_nonviol_sequence(neuron_count,:) = (nanmean(sdf.sequenceOnset(strcmp(event_table.cond_label,'nonviol'),:))-baseline_fr_mu_seq)./baseline_fr_std_seq;
        sdf_nonviol_reward(neuron_count,:) = (nanmean(sdf.reward(strcmp(event_table.cond_label,'nonviol'),:))-baseline_fr_mu_reward)./baseline_fr_std_reward;

    end
end



for i = 1:size(sdf_nonviol_sequence,1)
    sdf_nonviol_trialStart_smooth(i,:) = smooth(sdf_nonviol_trialStart(i,:),50);
    sdf_nonviol_sequence_smooth(i,:) = smooth(sdf_nonviol_sequence(i,:),50);
    sdf_nonviol_reward_smooth(i,:) = smooth(sdf_nonviol_reward(i,:),50);
end

%%



aud_mod_label = repmat({'Modulated'},length(aud_mod_neurons),1);
frontal_mod_label = repmat({'Modulated'},length(frontal_mod_neurons),1);
aud_nonmod_label = repmat({'Nonmodulated'},length(aud_nonmod_neurons),1);
frontal_nonmod_label = repmat({'Nonmodulated'},length(frontal_nonmod_neurons),1);


clear sequence_fig

sequence_fig(1,1)=gramm('x',ops.timewin,'y',sdf_nonviol_trialStart_smooth([aud_mod_neurons; aud_nonmod_neurons],:),'color',[aud_mod_label; aud_nonmod_label]);
sequence_fig(1,1).stat_summary();
sequence_fig(1,1).axe_property('XLim',[-100 500],'YLim',[-0.5 1]);
sequence_fig(1,1).geom_vline('xintercept',[0],'style','k-');
sequence_fig(1,1).no_legend;

sequence_fig(2,1)=gramm('x',ops.timewin,'y',sdf_nonviol_trialStart_smooth([frontal_mod_neurons; frontal_nonmod_neurons],:),'color',[frontal_mod_label; frontal_nonmod_label]);
sequence_fig(2,1).stat_summary();
sequence_fig(2,1).axe_property('XLim',[-100 500],'YLim',[-0.5 1]);
sequence_fig(2,1).geom_vline('xintercept',[0],'style','k-');
sequence_fig(2,1).no_legend;

sequence_fig(1,2)=gramm('x',ops.timewin,'y',sdf_nonviol_sequence_smooth([aud_mod_neurons; aud_nonmod_neurons],:),'color',[aud_mod_label; aud_nonmod_label]);
sequence_fig(1,2).stat_summary();
sequence_fig(1,2).axe_property('XLim',[-250 3000],'YLim',[-0.5 1]);
sequence_fig(1,2).geom_vline('xintercept',[0 563 1126 1689 2252],'style','k-');
sequence_fig(1,2).geom_vline('xintercept',[413 976 1539 2102 2665],'style','k--');
sequence_fig(1,2).no_legend;

sequence_fig(2,2)=gramm('x',ops.timewin,'y',sdf_nonviol_sequence_smooth([frontal_mod_neurons; frontal_nonmod_neurons],:),'color',[frontal_mod_label; frontal_nonmod_label]);
sequence_fig(2,2).stat_summary();
sequence_fig(2,2).axe_property('XLim',[-250 3000],'YLim',[-0.5 1]);
sequence_fig(2,2).geom_vline('xintercept',[0 563 1126 1689 2252],'style','k-');
sequence_fig(2,2).geom_vline('xintercept',[413 976 1539 2102 2665],'style','k--');
sequence_fig(2,2).no_legend;

sequence_fig(1,3)=gramm('x',ops.timewin,'y',sdf_nonviol_reward_smooth([aud_mod_neurons; aud_nonmod_neurons],:),'color',[aud_mod_label; aud_nonmod_label]);
sequence_fig(1,3).stat_summary();
sequence_fig(1,3).axe_property('XLim',[-250 1000],'YLim',[-0.5 2.5]);
sequence_fig(1,3).geom_vline('xintercept',[0],'style','k-');
sequence_fig(1,3).no_legend;

sequence_fig(2,3)=gramm('x',ops.timewin,'y',sdf_nonviol_reward_smooth([frontal_mod_neurons; frontal_nonmod_neurons],:),'color',[frontal_mod_label; frontal_nonmod_label]);
sequence_fig(2,3).stat_summary();
sequence_fig(2,3).axe_property('XLim',[-250 1000],'YLim',[-0.5 2.5]);
sequence_fig(2,3).geom_vline('xintercept',[0],'style','k-');
sequence_fig(2,3).no_legend;

% Figure setup    ////////////////////////////////////////////////////////
sequence_fig(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
sequence_fig(1,2).axe_property('XTick',[],'XColor',[1 1 1]);
sequence_fig(1,3).axe_property('XTick',[],'XColor',[1 1 1]);

% Figure layout ////////////////////////////////////////////////////////////
sequence_fig(1,1).set_layout_options...
    ('Position',[0.07 0.55 0.1 0.35],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

sequence_fig(2,1).set_layout_options...
    ('Position',[0.07 0.1 0.1 0.35],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

sequence_fig(1,2).set_layout_options...
    ('Position',[0.22 0.55 0.5 0.35],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

sequence_fig(2,2).set_layout_options...
    ('Position',[0.22 0.1 0.5 0.35],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

sequence_fig(1,3).set_layout_options...
    ('Position',[0.77 0.55 0.2 0.35],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

sequence_fig(2,3).set_layout_options...
    ('Position',[0.77 0.1 0.2 0.35],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

figure('Renderer', 'painters', 'Position', [100 100 1200 400]);
sequence_fig.draw




