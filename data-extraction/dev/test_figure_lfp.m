




for ch_i = 10
    clear test


    test(1,1)=gramm('x',ops.timewin,'y',lfp_aligned.(['lfp_' int2str(ch_i)]));
    test(1,1).stat_summary();
    test(1,1).axe_property('XLim',[-100 1000],'YLim',[-25 25]);


    test(1,1).geom_vline('xintercept',0,'style','k-');


    test(1,1).set_names('x','Time from event (ms)','y','Voltage');



    figure('Position',[100 100 1000 400]);
    test.draw();

end


