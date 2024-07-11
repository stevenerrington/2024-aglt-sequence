

%% IN DEV: 2024-07-08
% SE-  just developed some code to try and form an alignment point for the
% time of violation. Currently just single violation, and this is collapsed
% across all conditions

ops.timewin = -1000:5000;
ops.freq = [2 200];
session_i = 68;
datafile = ephysLog.session{session_i};
load(fullfile(dirs.mat_data,datafile))



aligntime = event_table.stimulusOnset_ms;
ops.timewin = -1000:5000;
ops.sdf_filter = 'Gauss';
[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);


sound_sdf_window = [-200:600];
sound_onset_ms = [0, 563, 1126, 1689, 2252];
zero_offset = abs(ops.timewin(1));

% single_unit_rastersdf_figure(event_table,spikes,ops)




%%
neuron_list = fieldnames(spikes.time);

for neuron_i = 1:length(neuron_list)
neuron_label = neuron_list{neuron_i};

sound_sdf = {};
count= 0;
for trial_i = 1:size(event_table,1)
    if ~strcmp(event_table.cond_label(trial_i),'error')
        for sound_i = 1:5
            count = count + 1;
            sound_sdf{count,1} = sdf.(neuron_label)(trial_i,zero_offset+sound_sdf_window+sound_onset_ms(sound_i));
            sound_sdf{count,2} = raster.(neuron_label){trial_i}-sound_onset_ms(sound_i);
            sound_sdf{count,3} = stimulusLog.(['sound_' int2str(sound_i) '_code']){event_table.cond_value(trial_i)};
            sound_sdf{count,4} = ['position_' int2str(sound_i)];
            sound_sdf{count,5} = event_table.cond_label(trial_i);
        end
    end
end

sound_list = {'A','C','D','F','G'};
sound_sdf_conc = struct();

sdf_in = cell2mat(sound_sdf(:,1));
raster_in = sound_sdf(:,2);
sound_cond = sound_sdf(:,3);
order_cond = sound_sdf(:,4);

xlim_vals = [-100 500];
ylim_vals = [0 10];
clear single_unit_fig

single_unit_fig(1,1)=gramm('x',raster_in);
single_unit_fig(1,1).geom_raster('geom',{'point'});
single_unit_fig(1,1).axe_property('XLim',xlim_vals);

single_unit_fig(2,1)=gramm('x',sound_sdf_window,'y',sdf_in);
single_unit_fig(2,1).stat_summary();
single_unit_fig(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);


single_unit_fig(1,1).axe_property('XTick',[],'XColor',[1 1 1]);

single_unit_fig(1,1).geom_vline('xintercept',0,'style','k-');
single_unit_fig(2,1).geom_vline('xintercept',0,'style','k-');


single_unit_fig(1,1).set_names('y','Trials');
single_unit_fig(2,1).set_names('x','Time from event (ms)','y','FR (spk/sec)');

single_unit_fig(1,1).set_layout_options...
    ('Position',[0.1 0.8 0.8 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(2,1).set_layout_options...
    ('Position',[0.1 0.1 0.8 0.6],... %Set the position in the figure (as in standard 'Position' axe property)
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);


single_unit_fig.set_title(neuron_label);
figure('Renderer', 'painters', 'Position', [100 100 700 600]);
single_unit_fig.draw();


end
















    figuren; 
    subplot(2,1,1); hold on
    for sound_list_i = 1:5
        sound_label = sound_list{sound_list_i};
        sound_sdf_conc.(['sound_' sound_label]) = cell2mat(sound_sdf(find(strcmp(sound_sdf(:,3),sound_label)),1));
        plot(sound_sdf_window, nanmean(sound_sdf_conc.(['sound_' sound_label])))
    end

    subplot(2,1,2); hold on
    plot(sound_sdf_window, nanmean(cell2mat(sound_sdf(:,1))))
    title(neuron_label)


%%
seq_cipher = {[3 7], [14], [1127];...
    [4 8], [15], [2253];...
    [1 5], [13], [1127];...
    [2 6], [16], [2253]};

neuron_list = fieldnames(spikes.time);


for neuron_i = 1:length(neuron_list)
    neuron_label = neuron_list{neuron_i};

    nonviol_sdf = []; viol_sdf = [];
    for seq_i = 1:size(seq_cipher,1)

        non_violation_trl = []; violation_trl = [];
        non_violation_trl = find(event_table.cond_value == seq_cipher{seq_i,1}(1) | event_table.cond_value == seq_cipher{seq_i,1}(2));
        violation_trl = find(event_table.cond_value == seq_cipher{seq_i,2}(1));

        analysis_window = [-100:500];
        zero_offset = abs(ops.timewin(1));

        nonviol_sdf = [nonviol_sdf; nanmean(sdf.(neuron_label)(non_violation_trl,zero_offset+analysis_window+seq_cipher{seq_i,3}(1)))];
        viol_sdf = [viol_sdf; nanmean(sdf.(neuron_label)(violation_trl,zero_offset+analysis_window+seq_cipher{seq_i,3}(1)))];
    end


    figuren; hold on
    plot(analysis_window,smooth(nanmean(nonviol_sdf),20),'g','LineWidth',2)
    plot(analysis_window,smooth(nanmean(viol_sdf),20),'r','LineWidth',2)
    title(neuron_label)
end



%% PATCH INTO EVENT TABLE FOR ALIGNMENT?
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

%%
aligntime = event_table.violation_ms;

ops.timewin = -1000:5000;

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);

names = fieldnames( spikes.time );

xlim_vals = [-100 300];
ylim_vals = [0 30];
clear test
% 49, 59 is good!


for ch_i = 1:length(names)

    ch = names{ch_i}(end-2:end);

    clear test

    test(1,1)=gramm('x',raster.(['DSP' ch]),'color',event_table.cond_label);
    test(1,1).geom_raster('geom',{'line'});
    test(1,1).axe_property('XLim',xlim_vals);


    test(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]),'color',event_table.cond_label);
    test(2,1).stat_summary();
    test(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);


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

end

