function single_unit_rastersdf_figure(spikes,ops)

aligntime = ops.aligntime;

ops.timewin = -1000:5000;
ops.sdf_filter = 'Gauss';

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);

names = ops.plot_ch;

xlim_vals = [-200 1200];
ylim_vals = [0 20];
close all
clear single_unit_fig


for ch_i = 1:length(names)

    ch = names{ch_i}(end-2:end);

    clear single_unit_fig

    single_unit_fig(1,1)=gramm('x',raster.(['DSP' ch]));
    single_unit_fig(1,1).geom_raster('geom',{'line'});
    single_unit_fig(1,1).axe_property('XLim',xlim_vals,'YLim',[0 size(raster.(['DSP' ch]),1)]);


    single_unit_fig(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]));
    single_unit_fig(2,1).stat_summary();
    single_unit_fig(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);


    single_unit_fig(3,1)=gramm('x',[-41:40],'y',spikes.waveform.(['WAV' ch]));
    single_unit_fig(3,1).stat_summary();
    single_unit_fig(3,1).axe_property('XLim',[-41 40],'YLim',[-35 35]);


    single_unit_fig(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
    single_unit_fig(3,1).axe_property('XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);


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

    single_unit_fig(3,1).set_layout_options...
        ('Position',[0.7 0.5 0.2 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
        'legend',false,...
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    single_unit_fig.set_title(['DSP' ch]);
    figure('Renderer', 'painters', 'Position', [100 100 700 600]);
    single_unit_fig.draw();

end
end