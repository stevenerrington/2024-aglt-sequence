function [encoding_flag, glm_beta, glm_sig, window_sdf, window_time] = glm_idpos_modulation(sound_info_in)

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

reg_tbl.firing_rate = win_fr;

%% Run GLM: trial type
% Identity model
pos_idx = strcmp(reg_tbl.condition,'nonviol') & reg_tbl.valid == 1 & strcmp(reg_tbl.sound,'C') & (strcmp(reg_tbl.order_pos,'position_2') | strcmp(reg_tbl.order_pos,'position_3') | strcmp(reg_tbl.order_pos,'position_4')  | strcmp(reg_tbl.order_pos,'position_5'));
id_idx = strcmp(reg_tbl.condition,'nonviol') &  reg_tbl.valid == 1 & strcmp(reg_tbl.order_pos,'position_3') & (strcmp(reg_tbl.sound,'C') | strcmp(reg_tbl.sound,'F') | strcmp(reg_tbl.sound,'G'));

reg_tbl.sound = categorical(reg_tbl.sound);
reg_tbl.order_pos = categorical(reg_tbl.order_pos);

reg_tbl_position = reg_tbl(pos_idx,:);
reg_tbl_identity = reg_tbl(id_idx,:);

u_t_mdl_identity = fitglm(reg_tbl_position, 'firing_rate ~ sound + exp_i', 'Distribution','Normal');
u_t_mdl_position = fitglm(reg_tbl_identity, 'firing_rate ~ order_pos + exp_i', 'Distribution','Normal');


end