function [glm_output, encoding_flag, encoding_beta, z_sdf] = violation_detect_glm(sdf_in, event_table_in)

%% Extract: get relevant data for GLM table
reg_tbl = table;

sound_info_in_sdf = sdf_in.sdf.violation(:,1000+[-200:800]);

reg_tbl.condition = event_table_in.event_table.cond_label;

mean_fr = nanmean(nanmean(sound_info_in_sdf(:,200+[-100:0])));
std_fr =  nanstd(nanmean(sound_info_in_sdf(:,200+[-100:0])));

%% Setup spike data into GLM
for trial_i = 1:size(sound_info_in_sdf,1)
    z_sdf(trial_i,:) = (sound_info_in_sdf(trial_i,:)-mean_fr)./std_fr;
    reg_tbl.exp_i(trial_i,1) = trial_i;
end

window_size = 100;
window_shift = 10;
[window_sdf, window_time] = movaverage_sdf(z_sdf, window_size, window_shift);

% 2023-08-26, 23h25: I tested this and it does the job correctly. The
% z-scored profile looks exactly like the raw profile, and the movaverage
% sdf looks like a smoothed version of the raw/zscored sdf.

[n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

baseline_win_idx = find(window_time >= -200 & window_time <= 0); % Find the relevant indicies for the timepoints of interest
analysis_win_idx = find(window_time >= 0 & window_time <= 600); % Find the relevant indicies for the timepoints of interest
win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

reg_tbl.win_fr = win_fr; % Add the firing rate over this whole window to the GLM table

cond_inc = [];
cond_inc = [1 5 14 2 6 15 3 7 13 4 8 16];

single_viol_trials = [];
single_viol_trials = find(ismember(event_table_in.event_table.cond_value,cond_inc) & ~isnan(event_table_in.event_table.rewardOnset_ms));

reg_tbl = reg_tbl(single_viol_trials,:);

%% Run GLM: trial type
glm_output.trial_type.sig_times = [];
glm_output.trial_type.beta_weights = [];
u_t_mdl = [];

% For each averaged time point
for timepoint_i = 1:n_times

    % Input the timepoint specific firing times
    reg_tbl.firing_rate = window_sdf(single_viol_trials,timepoint_i);

    % Convert 'condition' to a categorical variable if it's not already
    reg_tbl.condition = categorical(reg_tbl.condition);

    % Then include all these in your model
    u_t_mdl = fitlm(reg_tbl, 'firing_rate ~ exp_i + condition');

    % GLM output -------------------------------
    % - Violation
    glm_output.trial_type.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < .001; % trial type
    glm_output.trial_type.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.tStat(2); % trial type

end


%% Determine periods of significance

signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 50 ms of selectivity following SSD

% Trial-type
clear start sig_len
[start, sig_len, ~] = ZeroOnesCount(glm_output.trial_type.sig_times(:,analysis_win_idx)); % choice direction

encoding_flag = []; encoding_beta = [];

encoding_flag(1) = any(sig_len >= signal_detect_wins);

if encoding_flag(1) == 1
    sig_range = [];
    sig_range(1) = analysis_win_idx(1)+start(find(sig_len >= signal_detect_wins,1,'first'));
    sig_range(2) = sig_range(1)+sig_len(find(sig_len >= signal_detect_wins,1,'first'));


    encoding_beta(1) = nanmean(glm_output.trial_type.beta_weights(:,sig_range(1):sig_range(2)));
else
    encoding_beta(1) = NaN;
end

