filenames = getDirFilenames('C:\KIKUCHI-LOCAL\data\ephys\mat','troy');

for file_i = 1:length(filenames)
    try
        file_in = filenames{file_i};

        load(fullfile('C:\KIKUCHI-LOCAL\data\ephys\mat',file_in))

        % Patch faulty channel
        record_idx = find(strcmp(ephysLog.session,file_in(1:end-4)),1);
        fault_ch_idx = ephysLog.faulty_ch(record_idx);
        lfp = patch_fault_ch(lfp,fault_ch_idx);

        clear aligntime
        aligntime = event_table.stimulusOnset_ms;
        ops.timewin = -1000:5000;
        ops.freq = [2 30];
        lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

        channels = fieldnames(lfp_aligned);
        channels = channels(1:16);

        clear signal_out

        % - for each electrode
        for electrode_i = 1:length(channels)
            n_trials = size(lfp_aligned.(['lfp_' int2str(electrode_i)]),1);

            % - across each trial
            for trial_i = 1:n_trials

                % - get the LFP data for the given
                signal_in = lfp_aligned.(['lfp_' int2str(electrode_i)])(trial_i,:);

                % - save the relevant data in an output array for future use
                signal_out(electrode_i,:,trial_i) = signal_in; % nchans x trialtime x ntrials


            end
        end

        trial_average_lfp = [];
        trial_average_lfp = nanmean(signal_out,3);


        figuren('Renderer', 'painters', 'Position', [100 100 1200 400]);

        subplot(1,4,1); hold on
        for ch_i = 1:length(channels)
            if ismember(ch_i,[1:16])
                color_line = [204 0 102]./255;
            else
                color_line = [0 153 153]./255;
            end

            plot(ops.timewin,  trial_average_lfp(ch_i,:)+10*(ch_i-1),'color',color_line)
        end
        set(gca,'ydir', 'reverse')
        ylim([-15 (length(channels)*10)]); yticks([10*([1:32]-1)]); yticklabels(num2cell([1:32]))
        xlim([-100 500])
        title(file_in)


        subplot(1,4,2); hold on
        csd_out = D_CSD_BASIC(signal_out, 'spc', 0.2);
        ops.baseline = [-500:0];
        csd_out = adjustBaseline(csd_out,ops);
        data = nanmean(csd_out,3);

        % Plot CSD --------------------------------------------
        for ch_i = 1:16
            color_line = [204 0 102]./255;
            plot(ops.timewin,  data(ch_i,:)+20*(ch_i-1),'color',color_line)
        end
        set(gca,'ydir', 'reverse')
        ylim([-15 (16*20)]); yticks([20*([1:16]-1)]); yticklabels(num2cell([1:16]))
        xlim([-100 500])


        subplot(1,4,[3 4])
        csd_in = H_2DSMOOTH(data);

        % Get figure
        imagesc('XData',ops.timewin,'YData',1:size(csd_in,1),'CData',csd_in)
        xlim([-100 500])
        ylim([find(~isnan(csd_in(:,1)),1,'first') find(~isnan(csd_in(:,1)),1,'last')])
        set(gca,'YDir','Reverse')
        colorscale = flipud(cbrewer('div','RdBu',100));
        colorscale(colorscale<0) = 0;
        colormap(colorscale)
        colorbar
    end
end

