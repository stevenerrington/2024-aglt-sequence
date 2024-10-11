
sdf_color_list = cbrewer('qual','Dark2',5);

n_col = 5;
n_row = 5;
n_item_per_page = n_col*n_row;
count = 0;
page_count = 1;

f = figuren('Renderer', 'painters', 'Position', [100 100 1800 1000]);

for neuron_i = 1:size(spike_log,1)
    count = count + 1;
    
    if count > 25
        count = 1;
        page_count = page_count + 1;
        print(f,fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data-extraction\doc\sound_summary', ['sound_summary_pg' int2str(page_count) '.png']),'-dpng','-r1000');

        close all
        f = figuren('Renderer', 'painters', 'Position', [100 100 1800 1000]);

    end


    subplot(n_col,n_row,count); hold on
    sdf_in = [];
    sdf_in = cell2mat(sdf_soundAlign_data{neuron_i}(:,1));

    plot(ops.sound_sdf_window,smooth(nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'),:)),50),'color',sdf_color_list(1,:))
    plot(ops.sound_sdf_window,smooth(nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'),:)),50),'color',sdf_color_list(2,:))
    plot(ops.sound_sdf_window,smooth(nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'),:)),50),'color',sdf_color_list(3,:))
    plot(ops.sound_sdf_window,smooth(nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'),:)),50),'color',sdf_color_list(4,:))
    plot(ops.sound_sdf_window,smooth(nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'),:)),50),'color',sdf_color_list(5,:))
    plot(ops.sound_sdf_window,smooth(nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),:)),50),'color',[0 0 0 0.5])
    title(['Neuron ' int2str(neuron_i) '- glm sig: ' int2str(ismember(neuron_i,modulated_neurons))]); vline(0,'k-')



end