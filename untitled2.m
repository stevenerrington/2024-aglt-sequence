
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
ephysLog_all = import_exp_map();
ephysLog = clean_exp_map(ephysLog_all);

session_i = 1;

outfile_name = ephysLog.session{session_i}; % Processed file name
data_in = load(fullfile(dirs.mat_data,[outfile_name '.mat']));

aligntime_a = data_in.event_table.trialStart_ms ;
aligntime_b = data_in.event_table.stimulusOnset_ms;
aligntime_c = data_in.event_table.rewardOnset_ms;

ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';
[sdf_trialStart, raster_trialStart] = get_spikes_aligned(data_in.spikes,aligntime_a,ops);
[sdf_stimulus, raster_stimulus] = get_spikes_aligned(data_in.spikes,aligntime_b,ops);
[sdf_reward, raster_reward] = get_spikes_aligned(data_in.spikes,aligntime_c,ops);

trial_idx_stim_nonviol = find(strcmp(data_in.event_table.cond_label,'nonviol'));
trial_idx_reward = find(~isnan(aligntime_c));

session_allidx = find(strcmp(ephysLog.session{session_i},ephysLog_all.session));

sdf_session = SpkConvolver (round(data_in.spikes.time.(data_in.spk_info.unitDSP{spike_i})),...
    round(max(data_in.spikes.time.(data_in.spk_info.unitDSP{spike_i}))+5000), ops.sdf_filter);

for spike_i = 40:50

    figuren('Renderer', 'painters', 'Position', [100 100 1500 800]); hold on;

    if ismember(data_in.spk_info.site(spike_i),[1:16])
        area = ephysLog_all.area_label{session_allidx(1)};
    else ismember(data_in.spk_info.site(spike_i),[17:32])
        area = ephysLog_all.area_label{session_allidx(2)};
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
    set ( title_ax, 'visible', 'off')

    % Waveform
    nsubplot(5,4,[1 2],[2])
    plot(nanmean(data_in.spikes.waveform.(data_in.spk_info.unitWAV{spike_i})),'linewidth',2,'Color','k')
    set(gca,'xcolor','none','ycolor','none')

    nsubplot(5,4,[1 2],3)
    histogram(diff(data_in.spikes.time.(data_in.spk_info.unitDSP{spike_i})),1:1:50,'LineStyle','None','FaceColor',[14 80 102]./255)
    xlim([-5 50])
    xlabel('Time since last spike (ms)'); ylabel('Frequency')


    a = nsubplot(5,4,[ 4 5],[1]);
    plot(ops.timewin, smooth(nanmean(sdf_trialStart.(data_in.spk_info.unitDSP{spike_i})(trial_idx_reward,:)),50),'linewidth',1,'Color',[220 9 34]./255)
    xlim([-200 1000]); xlabel('Time from Trial Start (ms)'); ylabel('Firing rate (spk/sec)')
    vline(0,'k');

    b = nsubplot(5,4,[ 4 5],[2 3]);
    plot(ops.timewin, smooth(nanmean(sdf_stimulus.(data_in.spk_info.unitDSP{spike_i})(trial_idx_stim_nonviol,:)),50),'linewidth',1,'Color',[144 27 112]./255)
    xlim([-200 3500]); xlabel('Time from stimulus on (ms)'); ylabel('Firing rate (spk/sec)')
    vline(0,'k');

    c = nsubplot(5,4,[ 4 5],[4]);
    plot(ops.timewin, smooth(nanmean(sdf_reward.(data_in.spk_info.unitDSP{spike_i})(trial_idx_reward,:)),50),'linewidth',1,'Color',[41 47 86]./255)
    xlim([-200 1000]); xlabel('Time from reward (ms)'); ylabel('Firing rate (spk/sec)')
    vline(0,'k');

    linkaxes([a b c],'y')

end