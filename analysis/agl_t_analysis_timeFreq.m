%% Analysis
aligntime = events.stimulusOnset_ms;

% Setup LFP data for the time frequency analysis
ops.timewin = [-1000:7000];
ops.freq = [1 100];
lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

% Run time-frequency power extraction limited to trials of interest
ops.tf_trials = find(~isnan(events.rewardOnset_ms) & strcmp(events.cond_label,'viol'));
ops.ch_extract = [17:32];
tf_data = get_timefrequency(lfp_aligned, ops);

%% Figure

figuren('Renderer', 'painters', 'Position', [100 100 900 400]);
ax1 = nsubplot(1,1,1,1);
contourf(ops.timewin,tf_data.frequencies,nanmean(tf_data.power,3),40,'linecolor','none');
set(gca,'ytick',round(logspace(log10(tf_data.frequencies(1)),log10(tf_data.frequencies(end)),10)*100)/100,'yscale','log','xlim',[min(ops.timewin) max(ops.timewin)],'clim',[0 10])

colorscale = flipud(cbrewer('div','RdBu',100));
colorscale(colorscale<0) = 0;
colormap(colorscale)
colorbar
