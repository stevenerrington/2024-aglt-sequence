function [glm_output, encoding_flag_ord,...
    encoding_flag_rel, encoding_beta_ord, encoding_beta_rel] =...
    glm_position_relative(sound_info_in, transitional_probability)

%% Extract: get relevant data for GLM table
reg_tbl = table;

sound_info_in_sdf = cell2mat(sound_info_in(:,1));

reg_tbl.sound = string(sound_info_in(:,3));
reg_tbl.order_pos = string(sound_info_in(:,4));
reg_tbl.prev_element = string(sound_info_in(:,6));
reg_tbl.condition = string(sound_info_in(:,5));
reg_tbl.exp_i = [1:size(sound_info_in_sdf,1)]';





%% Setup spike data into GLM
mean_fr = nanmean(sound_info_in_sdf(:));
std_fr = nanstd(sound_info_in_sdf(:));

for trial_i = 1:size(sound_info_in_sdf,1)
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
analysis_win_idx = find(window_time >= 0 & window_time <= 600); % Find the relevant indicies for the timepoints of interest
win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

reg_tbl.win_fr = win_fr; % Add the firing rate over this whole window to the GLM table

%% Run GLM: trial type
clear glm_output
glm_output.trial_type.sig_times = [];
glm_output.trial_type.beta_weights = [];
u_t_mdl_position = [];
u_t_mdl_relative = [];

reg_tbl_trialtype = reg_tbl(strcmp(reg_tbl.condition,'nonviol') & strcmp(reg_tbl.sound,'C') ,:);

% For each averaged time point
for timepoint_i = 1:n_times

    % Input the timepoint specific firing times
    reg_tbl_trialtype.firing_rate = window_sdf(strcmp(reg_tbl.condition,'nonviol') & strcmp(reg_tbl.sound,'C'),timepoint_i);

    reg_tbl_trialtype.order_pos = categorical(reg_tbl_trialtype.order_pos);
    reg_tbl_trialtype.prev_element = categorical(reg_tbl_trialtype.prev_element);

    % Reorder categories to make "position_2" & "sound A" the reference
    reg_tbl_trialtype.order_pos = reordercats(reg_tbl_trialtype.order_pos, {'position_2', 'position_3', 'position_4', 'position_5'});
    reg_tbl_trialtype.prev_element = reordercats(reg_tbl_trialtype.prev_element, {'A', 'D', 'F', 'G'});

    % Then include all these in your model
    u_t_mdl_position = fitlm(reg_tbl_trialtype, 'firing_rate ~ exp_i + order_pos');
    u_t_mdl_relative = fitlm(reg_tbl_trialtype, 'firing_rate ~ exp_i + prev_element');


    % GLM output -------------------------------
    % - Ordinal Position
    glm_output.ordinal.sig_times(1,timepoint_i) = u_t_mdl_position.Coefficients.pValue(2) < .01; % Pos 3
    glm_output.ordinal.beta_weights(1,timepoint_i) = u_t_mdl_position.Coefficients.tStat(2); % Pos 3

    glm_output.ordinal.sig_times(2,timepoint_i) = u_t_mdl_position.Coefficients.pValue(3) < .01; % Pos 4
    glm_output.ordinal.beta_weights(2,timepoint_i) = u_t_mdl_position.Coefficients.tStat(3); % Pos 4

    glm_output.ordinal.sig_times(3,timepoint_i) = u_t_mdl_position.Coefficients.pValue(4) < .01; % Pos 5
    glm_output.ordinal.beta_weights(3,timepoint_i) = u_t_mdl_position.Coefficients.tStat(4); % Pos 5


    glm_output.ordinal.var_exp(1,timepoint_i) = u_t_mdl_position.Rsquared.Adjusted; % Pos 3
    glm_output.ordinal.LogLikelihood(1,timepoint_i) = u_t_mdl_position.LogLikelihood; % Pos 3

    % - Relative Position
    glm_output.relative.sig_times(1,timepoint_i) = u_t_mdl_relative.Coefficients.pValue(2) < .01; % Pos 3
    glm_output.relative.beta_weights(1,timepoint_i) = u_t_mdl_relative.Coefficients.tStat(2); % Pos 3

    glm_output.relative.sig_times(2,timepoint_i) = u_t_mdl_relative.Coefficients.pValue(3) < .01; % Pos 4
    glm_output.relative.beta_weights(2,timepoint_i) = u_t_mdl_relative.Coefficients.tStat(3); % Pos 4

    glm_output.relative.sig_times(3,timepoint_i) = u_t_mdl_relative.Coefficients.pValue(4) < .01; % Pos 5
    glm_output.relative.beta_weights(3,timepoint_i) = u_t_mdl_relative.Coefficients.tStat(4); % Pos 5

    glm_output.relative.var_exp(1,timepoint_i) = u_t_mdl_relative.Rsquared.Ordinary;
    glm_output.relative.LogLikelihood(1,timepoint_i) = u_t_mdl_relative.LogLikelihood;


end


%% Determine periods of significance

signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 50 ms of selectivity following SSD

% Trial-type
clear start_ord sig_len_ord start_rel sig_len_rel
for cond_i = 1:3
    [start_ord{cond_i}, sig_len_ord{cond_i}, ~] = ZeroOnesCount(glm_output.ordinal.sig_times(cond_i,analysis_win_idx)); % choice direction
    [start_rel{cond_i}, sig_len_rel{cond_i}, ~] = ZeroOnesCount(glm_output.relative.sig_times(cond_i,analysis_win_idx)); % choice direction
end

encoding_flag_ord = []; encoding_beta_ord = [];

for cond_i = 1:3
    encoding_flag_ord(1,cond_i) = any(sig_len_ord{cond_i} >= signal_detect_wins);

    if encoding_flag_ord(1,cond_i) == 1
        sig_range = [];
        sig_range(1) = analysis_win_idx(1)+start_ord{cond_i}(find(sig_len_ord{cond_i} >= signal_detect_wins,1,'first'));
        sig_range(2) = sig_range(1)+sig_len_ord{cond_i}(find(sig_len_ord{cond_i} >= signal_detect_wins,1,'first'));


        encoding_beta_ord(1,cond_i) = nanmean(glm_output.ordinal.beta_weights(cond_i,sig_range(1):sig_range(2)));
    else
        encoding_beta_ord(1,cond_i) = NaN;
    end
end


encoding_flag_rel = []; encoding_beta_rel = [];

for cond_i = 1:3
    encoding_flag_rel(1,cond_i) = any(sig_len_rel{cond_i} >= signal_detect_wins);

    if encoding_flag_rel(1,cond_i) == 1
        sig_range = [];
        sig_range(1) = analysis_win_idx(1)+start_rel{cond_i}(find(sig_len_rel{cond_i} >= signal_detect_wins,1,'first'));
        sig_range(2) = sig_range(1)+sig_len_rel{cond_i}(find(sig_len_rel{cond_i} >= signal_detect_wins,1,'first'));


        encoding_beta_rel(1,cond_i) = nanmean(glm_output.relative.beta_weights(cond_i,sig_range(1):sig_range(2)));
    else
        encoding_beta_rel(1,cond_i) = NaN;
    end
end
