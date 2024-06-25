clear all; clc
dirs = set_directories();
load(fullfile(dirs.mat_data,"troy-agl_t-2021-11-05.mat"))


% Improt data and get events
ops.event_port = 2;
events = get_agl_t_trials(event_table, ops);
aligntime = events.stimulusOnset_ms;

ops.timewin = [-1000:5000];
ops.freq = [1 100];
lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

%% 
ops.tf_trials = find(~isnan(events.rewardOnset_ms) & strcmp(events.cond_label,'nonviol'));
lfp_tf_out = get_timefrequency(lfp_aligned, ops);

%% Produce time frequency plot

mean_tf_sessions = nanmean(lfp_tf_out,3);

figuren;
ax1 = nsubplot(2,1,1,1);
contourf(lfp_tf.times,frequencies,mean_tf_sessions,25,'linecolor','none');
set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',[-1000 5000],'clim',[-10 10])

colorscale = flipud(cbrewer('div','RdBu',100));
colorscale(colorscale<0) = 0;
colormap(colorscale)
colorbar

ax2 = nsubplot(2,1,2,1);
theta_idx = find(frequencies > 1 & frequencies < 10);
plot(lfp_tf.times,nanmean(mean_tf_sessions(theta_idx,:)))