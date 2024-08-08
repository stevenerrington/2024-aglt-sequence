
%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- Kikuchi lab main script -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%} 

%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
[ephysLog_all, stimulusLog] = import_exp_map();
ephysLog = clean_exp_map(ephysLog_all);
load('session_audio_latency.mat')

tic
for session_i = 1 :size(ephysLog,1)

    outfile_name = ephysLog.session{session_i}; % Processed file name
    data_in = load(fullfile(dirs.mat_data,[outfile_name '.mat']));
    fprintf('Session %i of %i | %s \n', session_i, size(ephysLog,1), outfile_name)

    for trial_i = 1:size(data_in.event_table,1)
        data_in.event_table.stimulusOnset_ms(trial_i) = ...
            data_in.event_table.stimulusOnset_ms(trial_i) + ...
            session_audio_latency{session_i}(trial_i);
    end
    

    aligntime_a = data_in.event_table.trialStart_ms ;
    aligntime_b = data_in.event_table.stimulusOnset_ms;
    aligntime_c = data_in.event_table.rewardOnset_ms;

    ops.timewin = -1000:5000;
    ops.sdf_filter = 'PSP';
    [sdf_trialStart, raster_trialStart] = get_spikes_aligned(data_in.spikes,aligntime_a,ops);
    [sdf_stimulus, raster_stimulus] = get_spikes_aligned(data_in.spikes,aligntime_b,ops);
    [sdf_reward, raster_reward] = get_spikes_aligned(data_in.spikes,aligntime_c,ops);

    trial_idx_stim_nonviol = find(strcmp(data_in.event_table.cond_label,'nonviol'));
    trial_idx_stim_viol = find(strcmp(data_in.event_table.cond_label,'viol'));
    trial_idx_reward = find(~isnan(aligntime_c));
    trial_idx_stim_id = find(data_in.event_table.cond_value == 1);

    session_allidx = find(strcmp(ephysLog.session{session_i},ephysLog_all.session));




    onset_times = [ 0         563        1126        1689        2252];

    for spike_i = 1:length(fieldnames(data_in.spikes.time))

        sdf_session = SpkConvolver (round(data_in.spikes.time.(data_in.spk_info.unitDSP{spike_i})),...
            round(max(data_in.spikes.time.(data_in.spk_info.unitDSP{spike_i}))+5000), ops.sdf_filter);

        f = figuren('Renderer', 'painters', 'Position', [100 100 1500 800],'visible','off'); hold on;

        if ismember(data_in.spk_info.site(spike_i),[1:16]);
            area = ephysLog_all.area_label_recA{session_allidx(1)};
        else ismember(data_in.spk_info.site(spike_i),[17:32]);
            area = ephysLog_all.area_label_recA{session_allidx(2)};
        end

        % Session information
        title_ax = nsubplot(5,4,[1 2],[1]);
        text(-0.4,1.0,[data_in.spk_info.unitDSP{spike_i}],'FontSize',24,'FontWeight','bold');
        text(-0.4,0.8,['Area: ' area],'FontSize',15,'FontWeight','bold');
        text(-0.4,0.6,'Session Information','FontWeight','bold');
        text(-0.4,0.5,['Session: ' ephysLog.session{session_i}],'Interpreter', 'none');
        text(-0.4,0.4,['Monkey: ' ephysLog.monkey{session_i}],'Interpreter', 'none');
        text(-0.4,0.3,['Date: ' datestr(ephysLog.data(session_i))]);
        text(-0.4,0.2,['AP: ' int2str(ephysLog.nan_y(session_i))]);
        text(-0.4,0.1,['ML: ' int2str(ephysLog.nan_x(session_i))]);
        text(-0.4,0.0,['N Spikes: ' int2str(length(data_in.spikes.time.(data_in.spk_info.unitDSP{spike_i})))]);
        text(-0.4,-0.1,['Average FR: ' num2str(mean(sdf_session)) ' spk/sec']);
        set ( title_ax, 'visible', 'off');


        sdf_session = smooth(sdf_session, 50);

        % Waveform
        nsubplot(5,4,[1 2],[2]);
        plot(nanmean(data_in.spikes.waveform.(data_in.spk_info.unitWAV{spike_i})),'linewidth',2,'Color','k')
        set(gca,'xcolor','none','ycolor','none');

        nsubplot(5,4,[1 2],3);
        histogram(diff(data_in.spikes.time.(data_in.spk_info.unitDSP{spike_i})),1:1:50,'LineStyle','None','FaceColor',[14 80 102]./255)
        xlim([-5 50])
        xlabel('Time since last spike (ms)'); ylabel('Frequency')

        nsubplot(5,4,[1 2],4);
        time = []; time = 1:1:length(sdf_session);
        plot(movmean(time,10000), movmean(sdf_session,10000),'k')
        xlim([0 max(movmean(time,10000))])
        xlabel('Time from trial start (ms)'); ylabel('Firing rate (spk/sec)')

        a = nsubplot(5,4,[ 4 5],[1]);
        plot(ops.timewin, smooth(nanmean(sdf_trialStart.(data_in.spk_info.unitDSP{spike_i})(trial_idx_reward,:)),20),'linewidth',1,'Color',[220 9 34]./255)
        xlim([-200 1000]); xlabel('Time from Trial Start (ms)'); ylabel('Firing rate (spk/sec)')
        vline(0,'k');

        b = nsubplot(5,4,[ 4 5],[2 3]); hold on
        colorscale_blue = flipud(cbrewer('seq','Blues',8));
        colorscale_red = flipud(cbrewer('seq','Reds',8));

        colorscale_all = [colorscale_blue; colorscale_red];

        % for trial_type = 1:16
        %     trial_idx = []; trial_idx = find(data_in.event_table.cond_value == trial_type);
        %     plot(ops.timewin, smooth(nanmean(sdf_stimulus.(data_in.spk_info.unitDSP{spike_i})(trial_idx,:)),50),'linewidth',1,'Color',[colorscale_all(trial_type,:) 0.2])
        % end
        % plot(ops.timewin, smooth(nanmean(sdf_stimulus.(data_in.spk_info.unitDSP{spike_i})(trial_idx_stim_viol,:)),50),'linewidth',2,'Color',[colorscale_all(9,:) 1.0])

        plot(ops.timewin, smooth(nanmean(sdf_stimulus.(data_in.spk_info.unitDSP{spike_i})(trial_idx_reward,:)),20),'linewidth',2,'Color',[colorscale_all(1,:) 1.0])
        vline(onset_times,'r-')

        xlim([-200 3500]); xlabel('Time from stimulus on (ms)'); ylabel('Firing rate (spk/sec)')
        vline(0,'k');

        c = nsubplot(5,4,[ 4 5],[4]);
        plot(ops.timewin, smooth(nanmean(sdf_reward.(data_in.spk_info.unitDSP{spike_i})(trial_idx_reward,:)),20),'linewidth',1,'Color',[41 47 86]./255)
        xlim([-200 1000]); xlabel('Time from reward (ms)'); ylabel('Firing rate (spk/sec)')
        vline(0,'k');

        linkaxes([a b c],'y')
        % ylim([prctile(sdf_session, 25) prctile(sdf_session, 75)])


        figure_name = [outfile_name '_' data_in.spk_info.unitDSP{spike_i}];
        print(f,fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-laminar\data-extraction\doc\spk_summary_adjust', [figure_name '.png']),'-dpng','-r1000');
        close all
    end
end
toc