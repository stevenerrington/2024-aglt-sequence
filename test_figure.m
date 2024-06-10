
names = fieldnames( spikes.time );


for ch_i = 15:25

    ch = names{ch_i}(end-2:end);

    clear test

    test(1,1)=gramm('x',raster.(['DSP' ch]));
    test(1,1).geom_raster('geom',{'line'});
    test(1,1).axe_property('XLim',[-200 1000],'YLim',[0 size(raster.(['DSP' ch]),1)]);


    test(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]));
    test(2,1).stat_summary();
    test(2,1).axe_property('XLim',[-200 1000],'YLim',[0 10]);


    test(3,1)=gramm('x',[-41:40],'y',spikes.waveform.(['WAV' ch]));
    test(3,1).stat_summary();
    test(3,1).axe_property('XLim',[-10 20],'YLim',[-25 25]);


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
        'legend',false,...
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    test(3,1).set_layout_options...
        ('Position',[0.7 0.5 0.2 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
        'legend',false,...
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    figure('Position',[100 100 700 600]);
    test.draw();

end