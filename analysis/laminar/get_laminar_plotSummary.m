%% Set figure parameters
fprintf('Generating figure   \n')

% Colormap
clear color_heatmap
color_heatmap = flipud(cbrewer('div','RdBu',100));
color_heatmap(color_heatmap<0) = 0;

%% Create figure window
figuren('Renderer', 'painters', 'Position', [100 100 1800 600]); hold on;

area_labels = fieldnames(laminar_info);

for area_i = 1:2
    area_name = area_labels{area_i};

    switch area_name
        case 'auditory'
            channels = 1:16;
        case 'vlpfc'
            channels = [17:22,24:32];
    end

    %% Subplots

    % LFP -----------------------------------------------------------------
    trial_average_lfp = nanmean(laminar_info.(area_name).lfp,3);

    ax_lfp = nsubplot(2,4,area_i,1);
    for ch_i = 1:length(channels)
        plot(ops.timewin,  trial_average_lfp(ch_i,:)+10*(ch_i-1),'k')
    end
    set(gca,'ydir', 'reverse')
    ylim([-15 (length(channels)*10)]); yticks([])
    xlim([-1000 5000])
    ylabel(area_name,'fontsize', 12, 'fontweight', 'bold')

    % CSD -----------------------------------------------------------------
    ax_csd = nsubplot(2,4,area_i,2);
    plot_csd(nanmean(laminar_info.(area_name).csd,3), ops.timewin, [-100 500], ax_csd)
    set(gca,'ydir', 'reverse')
    cb=colorbar; ylabel(cb, 'nA/mm3');
    xlabel('Time (ms)');ylabel('Depth');

    % Spectrolaminar power ------------------------------------------------
    ax_specpower = nsubplot(2,4,area_i,3);
    imagesc(laminar_info.(area_name).power.normalized);set(gca, 'YDir', 'reverse');
    xlim([1 size(laminar_info.(area_name).power.normalized,2)]); ylim([1 size(laminar_info.(area_name).power.normalized,1)]);
    xlabel('Frequency (Hz)');ylabel('Channel Number');
    colormap(color_heatmap)
    hline(laminar_info.(area_name).flip.crossoverchannel,'k')
    cb=colorbar; ylabel(cb, 'Relative Power');

    % Spectrolaminar summary ----------------------------------------------
    nsubplot(2,4,area_i,4);
    plot(mean(laminar_info.(area_name).power.normalized(:,10:19),2), 1:length(channels), 'b', 'LineWidth', 2); % alpha/beta
    plot(mean(laminar_info.(area_name).power.normalized(:,75:150),2), 1:length(channels), 'r', 'LineWidth', 2); % gamma
    set(gca, 'ydir', 'reverse')
    xlim([0 1]); ylim([1 length(channels)]); xlabel('Relative Power'); ylabel('Channel Number');
    hline(laminar_info.(area_name).flip.crossoverchannel,'k')

    % Other -----------------------------------------------------------------
end


    sgtitle( outfile_name , 'interpreter', 'none')
