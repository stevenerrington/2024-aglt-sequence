aligntime = event_table.rewardOnset_ms;

ops.timewin = -5000:5000;
ops.sdf_filter = 'PSP';

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);


clear pop_av_test

for ch_i = 1:length(names)
    ch = names{ch_i}(end-2:end);

    bl_mean = nanmean(nanmean(sdf.(['DSP' ch])(:,[-500:0]+1000)));
    bl_sd = nanstd(nanmean(sdf.(['DSP' ch])(:,[-500:0]+1000)));
    sdf_in = (sdf.(['DSP' ch]) - bl_mean)./bl_sd;

    test_sdf = [];



    %figuren('Renderer', 'painters', 'Position', [100 100 700 600]); hold on
    %plot(ops.timewin,nanmean(sdf_in),'color',[0 0 1],'LineWidth',1.5)

    pop_av_test(ch_i,:) = nanmean(sdf_in);

end



figuren('Renderer', 'painters', 'Position', [100 100 700 600]); hold on
subplot(2,1,1)
plot(ops.timewin,nanmean(pop_av_test(1:20,:)),'color',[0 0 1],'LineWidth',1.5)
subplot(2,1,2)
plot(ops.timewin,nanmean(pop_av_test(21:end,:)),'color',[0 0 1],'LineWidth',1.5)
