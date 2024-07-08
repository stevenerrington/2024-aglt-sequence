

%% IN DEV: 2024-07-08
% SE-  just developed some code to try and form an alignment point for the
% time of violation. Currently just single violation, and this is collapsed
% across all conditions

ops.timewin = -1000:5000;
ops.freq = [2 200];
session_i = 23;
datafile = ephysLog.session{session_i};
load(fullfile(dirs.mat_data,datafile))

aligntime = event_table.stimulusOnset_ms;
ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';
[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);

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
        viol_sdf = [viol_sdf; nanmean(sdf.(neuron_label)(violation_trl,zero_offset+analysis_window+seq_cipher{seq_i,3}(1))];
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
ops.sdf_filter = 'Gauss';

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);

names = fieldnames( spikes.time );

xlim_vals = [-100 500];
ylim_vals = [0 15];
close all
clear test
% 49, 59 is good!


for ch_i = 49

    ch = names{ch_i}(end-2:end);

    clear test

    test(1,1)=gramm('x',raster.(['DSP' ch]),'color',event_table.cond_label);
    test(1,1).geom_raster('geom',{'line'});
    test(1,1).axe_property('XLim',xlim_vals);


    test(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]),'color',event_table.cond_label);
    test(2,1).stat_summary();
    test(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);


    test(3,1)=gramm('x',[-41:40],'y',spikes.waveform.(['WAV' ch]));
    test(3,1).stat_summary();
    test(3,1).axe_property('XLim',[-41 40],'YLim',[-35 35]);


    test(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
    test(3,1).axe_property('XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);


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

    test(3,1).set_layout_options...
        ('Position',[0.7 0.5 0.2 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
        'legend',false,...
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    test.set_title(['DSP' ch]);
    figure('Renderer', 'painters', 'Position', [100 100 700 600]);
    test.draw();

end

