




for ch_i = 1:32

    clear test


    test(1,1)=gramm('x',ops.timewin,'y',lfp.(['lfp_' int2str(ch_i)]));
    test(1,1).stat_summary();
    test(1,1).axe_property('XLim',[-200 1000],'YLim',[-20 20]);


    test(1,1).geom_vline('xintercept',0,'style','k-');


    test(1,1).set_names('x','Time from event (ms)','y','Voltage');



    figure('Position',[100 100 700 600]);
    test.draw();

end