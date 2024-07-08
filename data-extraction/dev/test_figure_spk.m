aligntime = event_table.violation_ms;

ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);




names = fieldnames( spikes.time );

xlim_vals = [-100 350];
ylim_vals = [0 25];
close all
clear test
for ch_i = 66

    ch = names{ch_i}(end-2:end);

    clear test

    test(1,1)=gramm('x',raster.(['DSP' ch]),'color',event_table.cond_label);
    test(1,1).geom_raster('geom',{'line'});
    test(1,1).axe_property('XLim',xlim_vals,'YLim',[0 size(raster.(['DSP' ch]),1)]);


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






%%










% Get list of neurons
names = fieldnames( spikes.time );

% Define times for figures
xlim_vals = [-1000 5000];
ylim_vals = [0 30];

for ch_i = 1:length(names)
    clear spk_figure_idv fig

    % Get channel label
    ch = names{ch_i}(4:end);

    % Plot raster
    spk_figure_idv(1,1)=gramm('x',raster.(['DSP' ch]),'color',events.cond_label);
    spk_figure_idv(1,1).geom_raster('geom',{'line'});
    spk_figure_idv(1,1).axe_property('XLim',xlim_vals,'YLim',[0 size(raster.(['DSP' ch]),1)]);

    % Plot SDF
    spk_figure_idv(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]),'color',events.cond_label);
    spk_figure_idv(2,1).stat_summary();
    spk_figure_idv(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);

    % Plot waveform average
    spk_figure_idv(3,1)=gramm('x',[1:length(spikes.waveform.(['WAV' ch]))],'y',spikes.waveform.(['WAV' ch]));
    spk_figure_idv(3,1).stat_summary();
    spk_figure_idv(3,1).axe_property('XLim',[-41 40],'YLim',[-35 35]);

    % Figure properties
    spk_figure_idv(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
    spk_figure_idv(3,1).axe_property('XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);
    spk_figure_idv(1,1).geom_vline('xintercept',0,'style','k-');
    spk_figure_idv(2,1).geom_vline('xintercept',0,'style','k-');
    spk_figure_idv(1,1).set_names('y','Trials');
    spk_figure_idv(2,1).set_names('x','Time from event (ms)','y','FR (spk/sec)');

    % Figure layout
    spk_figure_idv(1,1).set_layout_options...
        ('Position',[0.1 0.8 0.8 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
        'legend',false,...
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    spk_figure_idv(2,1).set_layout_options...
        ('Position',[0.1 0.1 0.8 0.6],... %Set the position in the figure (as in standard 'Position' axe property)
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    spk_figure_idv(3,1).set_layout_options...
        ('Position',[0.7 0.5 0.2 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
        'legend',false,...
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    spk_figure_idv.set_title(['DSP' ch]);

    % Draw figure
    fig = figure('Renderer', 'painters', 'Position', [100 100 700 600]);
    spk_figure_idv.draw();

    % Save figure
    set(fig,'renderer','painters','Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig,['C:\KIKUCHI-LOCAL\script\kikuchi-data\data-extraction\doc\troy-agl_t-2021-11-05-Kilosort-',['DSP' ch],'.pdf'],'-r400','-bestfit','-dpdf')
    close all
end
