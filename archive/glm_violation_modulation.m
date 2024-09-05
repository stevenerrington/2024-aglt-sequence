function [glm_output, encoding_flag, encoding_beta, z_sdf] = glm_violation_modulation(sound_info_in)

%% Extract: get relevant data for GLM table
reg_tbl = table;

sound_info_in_sdf = cell2mat(sound_info_in(:,1));

reg_tbl.sound = string(sound_info_in(:,3));
reg_tbl.order_pos = string(sound_info_in(:,4));
reg_tbl.condition = string(sound_info_in(:,5));
reg_tbl.exp_i = [1:size(sound_info_in_sdf,1)]';


%% Setup spike data into GLM
for trial_i = 1:size(sound_info_in_sdf,1)
    mean_fr = nanmean(nanmean(sound_info_in_sdf(strcmp(reg_tbl.sound,reg_tbl.sound{trial_i}),200+[-100:0])));
    std_fr = std(nanmean(sound_info_in_sdf(strcmp(reg_tbl.sound,reg_tbl.sound{trial_i}),200+[-100:0])));

    z_sdf(trial_i,:) = (sound_info_in_sdf(trial_i,:)-mean_fr)./std_fr;
end

%%
window_size = 50;
window_shift = 10;
[window_sdf, window_time] = movaverage_sdf(z_sdf, window_size, window_shift);
% 2023-08-26, 23h25: I tested this and it does the job correctly. The
% z-scored profile looks exactly like the raw profile, and the movaverage
% sdf looks like a smoothed version of the raw/zscored sdf.

[n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

baseline_win_idx = find(window_time >= -200 & window_time <= 0); % Find the relevant indicies for the timepoints of interest
analysis_win_idx = find(window_time >= 0 & window_time <= 400); % Find the relevant indicies for the timepoints of interest
win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

reg_tbl.win_fr = win_fr; % Add the firing rate over this whole window to the GLM table

%% Run GLM: trial type
glm_output.violation.sig_times = [];
glm_output.violation.beta_weights = [];
u_t_mdl_all = []; u_t_mdl_c = []; u_t_mdl_f = []; u_t_mdl_g = [];

reg_tbl_trialtype_all = reg_tbl(strcmp(reg_tbl.sound,'C') | strcmp(reg_tbl.sound,'F') | strcmp(reg_tbl.sound,'G') ,:);
reg_tbl_trialtype_c = reg_tbl(strcmp(reg_tbl.sound,'C') ,:);
reg_tbl_trialtype_f = reg_tbl(strcmp(reg_tbl.sound,'F') ,:);
reg_tbl_trialtype_g = reg_tbl(strcmp(reg_tbl.sound,'G') ,:);


% For each averaged time point
for timepoint_i = 1:n_times

    % Input the timepoint specific firing times
    reg_tbl_trialtype_all.firing_rate = window_sdf(strcmp(reg_tbl.sound,'C') | strcmp(reg_tbl.sound,'F') | strcmp(reg_tbl.sound,'G'),timepoint_i);
    reg_tbl_trialtype_c.firing_rate = window_sdf(strcmp(reg_tbl.sound,'C'),timepoint_i);
    reg_tbl_trialtype_f.firing_rate = window_sdf(strcmp(reg_tbl.sound,'F'),timepoint_i);
    reg_tbl_trialtype_g.firing_rate = window_sdf(strcmp(reg_tbl.sound,'G'),timepoint_i);

    % Convert 'sound' to a categorical variable if it's not already
    reg_tbl_trialtype_all.condition = categorical(reg_tbl_trialtype_all.condition);
    reg_tbl_trialtype_c.condition = categorical(reg_tbl_trialtype_c.condition);
    reg_tbl_trialtype_f.condition = categorical(reg_tbl_trialtype_f.condition);
    reg_tbl_trialtype_g.condition = categorical(reg_tbl_trialtype_g.condition);

    % Reorder categories to make "G" the reference
    reg_tbl_trialtype_all.condition = reordercats(reg_tbl_trialtype_all.condition, {'nonviol', 'viol'});
    reg_tbl_trialtype_c.condition = reordercats(reg_tbl_trialtype_c.condition, {'nonviol', 'viol'});
    reg_tbl_trialtype_f.condition = reordercats(reg_tbl_trialtype_f.condition, {'nonviol', 'viol'});
    reg_tbl_trialtype_g.condition = reordercats(reg_tbl_trialtype_g.condition, {'nonviol', 'viol'});

    % Then include all these in your model
    % u_t_mdl_all = fitlm(reg_tbl_trialtype_all, 'firing_rate ~ exp_i + condition');
    u_t_mdl_c = fitlm(reg_tbl_trialtype_c, 'firing_rate ~ exp_i + condition');
    u_t_mdl_f = fitlm(reg_tbl_trialtype_f, 'firing_rate ~ exp_i + condition');
    u_t_mdl_g = fitlm(reg_tbl_trialtype_g, 'firing_rate ~ exp_i + condition');

    % GLM output -------------------------------
    % - Sound A
    glm_output.violation.sig_times(1,timepoint_i) = u_t_mdl_c.Coefficients.pValue(2) < .01; % trial type
    glm_output.violation.sig_times(2,timepoint_i) = u_t_mdl_f.Coefficients.pValue(2) < .01; % trial type
    glm_output.violation.sig_times(3,timepoint_i) = u_t_mdl_g.Coefficients.pValue(2) < .01; % trial type

    glm_output.violation.beta_weights(1,timepoint_i) = u_t_mdl_c.Coefficients.tStat(2); % trial type
    glm_output.violation.beta_weights(2,timepoint_i) = u_t_mdl_f.Coefficients.tStat(2); % trial type
    glm_output.violation.beta_weights(3,timepoint_i) = u_t_mdl_g.Coefficients.tStat(2); % trial type

    %glm_output.violation_all.sig_times(1,timepoint_i) = u_t_mdl_all.Coefficients.pValue(2) < .01; % trial type
    %glm_output.violation_all.beta_weights(1,timepoint_i) = u_t_mdl_all.Coefficients.tStat(2); % trial type

end


%% Determine periods of significance

signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 50 ms of selectivity following SSD

% Trial-type
clear sig_len
[~, sig_len, ~] = ZeroOnesCount(sum(glm_output.violation.sig_times(:,analysis_win_idx)) > 0); % choice direction


encoding_flag = []; encoding_beta = [];

encoding_flag(1,1) = any(sig_len >= signal_detect_wins);
encoding_beta(1,:) = nanmean(glm_output.violation.beta_weights(:,:));




end