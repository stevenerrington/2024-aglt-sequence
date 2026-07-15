function [encoding_flag] = ttest_violation_modulation(sound_info_in)

%% Extract: get relevant data for GLM table
reg_tbl = table;

sound_info_in_sdf = cell2mat(sound_info_in(:,1));

reg_tbl.sound = string(sound_info_in(:,3));
reg_tbl.order_pos = string(sound_info_in(:,4));
reg_tbl.condition = string(sound_info_in(:,5));
reg_tbl.exp_i = [1:size(sound_info_in_sdf,1)]';
reg_tbl.valid = cell2mat(sound_info_in(:,9));


%% Setup spike data into GLM
for trial_i = 1:size(sound_info_in_sdf,1)
    mean_fr = nanmean(nanmean(sound_info_in_sdf(strcmp(reg_tbl.sound,reg_tbl.sound{trial_i}),200+[-100:0])));
    std_fr = std(nanmean(sound_info_in_sdf(strcmp(reg_tbl.sound,reg_tbl.sound{trial_i}),200+[-100:0])));

    z_sdf(trial_i,:) = (sound_info_in_sdf(trial_i,:)-mean_fr)./std_fr;
end

window_size = 100;
window_shift = 10;
[window_sdf, window_time] = movaverage_sdf(z_sdf, window_size, window_shift);

% 2023-08-26, 23h25: I tested this and it does the job correctly. The
% z-scored profile looks exactly like the raw profile, and the movaverage
% sdf looks like a smoothed version of the raw/zscored sdf.

[n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

baseline_win_idx = find(window_time >= -200 & window_time <= 0); % Find the relevant indicies for the timepoints of interest
analysis_win_idx = find(window_time >= 0 & window_time <= 413); % Find the relevant indicies for the timepoints of interest
win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

reg_tbl.win_fr = win_fr; % Add the firing rate over this whole window to the GLM table

%% Run GLM: trial type
glm_output.trial_type.sig_times = [];
glm_output.trial_type.beta_weights = [];
u_t_mdl = [];

valid_trial_idx = [];
valid_trial_idx = find(reg_tbl.valid == 1 );

reg_tbl_trialtype = [];
reg_tbl_trialtype = reg_tbl(valid_trial_idx ,:);

sig_flag = []; p_value = [];

% For each averaged time point
for timepoint_i = 1:n_times

    % Input the timepoint specific firing times
    reg_tbl_trialtype.firing_rate = window_sdf(valid_trial_idx,timepoint_i);

    alpha_level = 0.01;

    [sig_flag(1,timepoint_i), p_value(1,timepoint_i)] =...
        ttest2(reg_tbl_trialtype.firing_rate(strcmp(reg_tbl_trialtype.condition,'nonviol'),:),...
        reg_tbl_trialtype.firing_rate(strcmp(reg_tbl_trialtype.condition,'viol'),:),...
        'Alpha',alpha_level,'Tail','both');

end


%% Determine periods of significance

signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 50 ms of selectivity following element

% Trial-type
clear start sig_len
[start, sig_len, ~] = ZeroOnesCount(sig_flag(1,analysis_win_idx));

encoding_flag = [];

encoding_flag = any(sig_len >= signal_detect_wins);


end