function [glm_output, encoding_flag, encoding_beta, z_sdf] = glm_sound_modulation(sound_info_in)

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

% For each averaged time point
for timepoint_i = 1:n_times

    % Input the timepoint specific firing times
    reg_tbl_trialtype.firing_rate = window_sdf(valid_trial_idx,timepoint_i);
    % Convert 'sound' to a categorical variable if it's not already
    reg_tbl_trialtype.sound = categorical(reg_tbl_trialtype.sound);

    % Reorder categories to make "G" the reference
    reg_tbl_trialtype.sound = reordercats(reg_tbl_trialtype.sound, {'Baseline', 'A', 'C', 'D', 'F','G'});

    % Then include all these in your model
    u_t_mdl = fitlm(reg_tbl_trialtype, 'firing_rate ~ sound');

    % GLM output -------------------------------
    % - Sound A
    glm_output.trial_type.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < (.0001/5); % trial type
    glm_output.trial_type.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.tStat(2); % trial type
    
    % - Sound C
    glm_output.trial_type.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < (.0001/5); % trial type
    glm_output.trial_type.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % trial type
    
    % - Sound D
    glm_output.trial_type.sig_times(3,timepoint_i) = u_t_mdl.Coefficients.pValue(4) < (.0001/5); % trial type
    glm_output.trial_type.beta_weights(3,timepoint_i) = u_t_mdl.Coefficients.tStat(4); % trial type
    
    % - Sound F
    glm_output.trial_type.sig_times(4,timepoint_i) = u_t_mdl.Coefficients.pValue(5) < (.0001/5); % trial type
    glm_output.trial_type.beta_weights(4,timepoint_i) = u_t_mdl.Coefficients.tStat(5); % trial type

    % - Sound G
    glm_output.trial_type.sig_times(5,timepoint_i) = u_t_mdl.Coefficients.pValue(6) < (.0001/5); % trial type
    glm_output.trial_type.beta_weights(5,timepoint_i) = u_t_mdl.Coefficients.tStat(6); % trial type
end


%% Determine periods of significance

signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 50 ms of selectivity following SSD

% Trial-type
clear start sig_len
for cond_i = 1:5
    [start{cond_i}, sig_len{cond_i}, ~] = ZeroOnesCount(glm_output.trial_type.sig_times(cond_i,analysis_win_idx)); % choice direction
end

encoding_flag = []; encoding_beta = [];

for cond_i = 1:5
    encoding_flag(1,cond_i) = any(sig_len{cond_i} >= signal_detect_wins);

    if encoding_flag(1,cond_i) == 1
        sig_range = [];
        sig_range(1) = analysis_win_idx(1)+start{cond_i}(find(sig_len{cond_i} >= signal_detect_wins,1,'first'));
        sig_range(2) = sig_range(1)+sig_len{cond_i}(find(sig_len{cond_i} >= signal_detect_wins,1,'first'));


        encoding_beta(1,cond_i) = nanmean(glm_output.trial_type.beta_weights(cond_i,sig_range(1):sig_range(2)));
    else
        encoding_beta(1,cond_i) = NaN;
    end
end


end