function [encoding_flag, glm_beta, glm_sig, window_sdf, window_time] = glm_sound_modulation(sound_info_in)

%% Extract: get relevant data for GLM table
reg_tbl = table;

sound_info_in_sdf = cell2mat(sound_info_in(:,1));

reg_tbl.sound = string(sound_info_in(:,3));
reg_tbl.order_pos = string(sound_info_in(:,4));
reg_tbl.condition = string(sound_info_in(:,5));
reg_tbl.exp_i = [1:size(sound_info_in_sdf,1)]';
reg_tbl.valid = cell2mat(sound_info_in(:,9));

%% Setup spike data into GLM
baseline_idx = 200+(-100:0);
mean_fr = nanmean(nanmean(sound_info_in_sdf(:,baseline_idx),2));
std_fr = nanstd(nanmean(sound_info_in_sdf(:,baseline_idx),2));
z_sdf = (sound_info_in_sdf - mean_fr) ./ std_fr;

window_size = 100;
window_shift = 10;
[window_sdf, window_time] = movaverage_sdf(sound_info_in_sdf, window_size, window_shift);
% 2023-08-26, 23h25: I tested this and it does the job correctly. The
% z-scored profile looks exactly like the raw profile, and the movaverage
% sdf looks like a smoothed version of the raw/zscored sdf.

[n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

baseline_win_idx = find(window_time >= -200 & window_time <= 0); % Find the relevant indicies for the timepoints of interest
analysis_win_idx = find(window_time >= 0 & window_time <= 413); % Find the relevant indicies for the timepoints of interest
win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

reg_tbl.win_fr = win_fr; % Add the firing rate over this whole window to the GLM table


%% Run GLM: trial type
p_cutoff = 0.05; % for FDR, use 0.05 rather than 0.001

glm_beta = nan(10, n_times);
glm_pval = nan(10, n_times);

valid_trial_idx = find(reg_tbl.valid == 1);
reg_tbl_trialtype = reg_tbl(valid_trial_idx ,:);

for timepoint_i = 1:n_times
    reg_tbl_trialtype.firing_rate = window_sdf(valid_trial_idx,timepoint_i);
    reg_tbl_trialtype.sound = categorical(reg_tbl_trialtype.sound);
    reg_tbl_trialtype.order_pos = categorical(reg_tbl_trialtype.order_pos);

    reg_tbl_trialtype.sound = reordercats(reg_tbl_trialtype.sound, {'Baseline', 'A', 'C', 'D', 'F','G'});

    % Identity model
    u_t_mdl_identity = fitglm(reg_tbl_trialtype, 'firing_rate ~ sound + exp_i', 'Distribution','Normal');
    % Position model
    u_t_mdl_position = fitglm(reg_tbl_trialtype, 'firing_rate ~ order_pos + exp_i', 'Distribution','Normal');

    % Extract betas and p-values
    glm_beta(1:5, timepoint_i) = u_t_mdl_identity.Coefficients.Estimate(2:6);
    glm_pval(1:5, timepoint_i) = u_t_mdl_identity.Coefficients.pValue(2:6);

    glm_beta(6:10, timepoint_i) = u_t_mdl_position.Coefficients.Estimate(2:6);
    glm_pval(6:10, timepoint_i) = u_t_mdl_position.Coefficients.pValue(2:6);
end

%% --- Apply FDR correction across time for each regressor ---
glm_sig = false(size(glm_pval));

for cond_i = 1:10
    this_p = glm_pval(cond_i, :);
    this_q = mafdr(this_p, 'BHFDR', true); % Benjamini-Hochberg correction
    glm_sig(cond_i, :) = this_q < p_cutoff;
end

%% Determine periods of significance

signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 50 ms of selectivity following element

% Trial-type
clear start sig_len
for cond_i = 1:10
    [start{cond_i}, sig_len{cond_i}, ~] = ZeroOnesCount(glm_sig(cond_i,analysis_win_idx));
end

encoding_flag = []; encoding_beta = [];

for cond_i = 1:10
    encoding_flag(1,cond_i) = any(sig_len{cond_i} >= signal_detect_wins);
end


end