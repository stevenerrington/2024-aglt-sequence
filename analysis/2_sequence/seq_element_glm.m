function [glm_sequence_out, anova_sequence_out] = seq_element_glm(sound_info_in, transitional_probability)

%% Extract: get relevant data for GLM table
reg_tbl = table;

sound_info_in_sdf = cell2mat(sound_info_in(:,1));

reg_tbl.sound = string(sound_info_in(:,3));
reg_tbl.order_pos = string(sound_info_in(:,4));
reg_tbl.condition = string(sound_info_in(:,5));
reg_tbl.exp_i = [1:size(sound_info_in_sdf,1)]';


for trial_i = 1:size(reg_tbl)

    if trial_i == 1 | strcmp(reg_tbl.sound(trial_i),'Baseline')
        reg_tbl.backward_prob(trial_i) = NaN;
        reg_tbl.forward_prob(trial_i) = NaN;
    else
        curr_trial_ele = find(strcmp(transitional_probability.elements, reg_tbl.sound(trial_i)));
        prev_trial_ele = find(strcmp(transitional_probability.elements, reg_tbl.sound(trial_i-1)));

        try
            reg_tbl.backward_prob(trial_i) = transitional_probability.backward_prob(prev_trial_ele,curr_trial_ele);
            reg_tbl.forward_prob(trial_i) = transitional_probability.forward_prob(prev_trial_ele,curr_trial_ele);
        catch
            reg_tbl.backward_prob(trial_i) = NaN;
            reg_tbl.forward_prob(trial_i) = NaN;
        end

    end

end

%% Setup spike data into GLM
mean_fr = nanmean(sound_info_in_sdf(:));
std_fr = nanstd(sound_info_in_sdf(:));

for trial_i = 1:size(sound_info_in_sdf,1)
    z_sdf(trial_i,:) = (sound_info_in_sdf(trial_i,:)-mean_fr)./std_fr;
end

window_size = 100;
window_shift = 10;
[window_sdf, window_time] = movaverage_sdf(z_sdf, window_size, window_shift);

[n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

baseline_win_idx = find(window_time >= -200 & window_time <= 0); % Find the relevant indicies for the timepoints of interest
analysis_win_idx = find(window_time >= 0 & window_time <= 500); % Find the relevant indicies for the timepoints of interest
win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

reg_tbl.firing_rate = win_fr; % Add the firing rate over this whole window to the GLM table

%% Curate tables for specific trials and conditions

% Initialise different model tables/data
reg_tbl_all = reg_tbl(strcmp(reg_tbl.condition,'nonviol') ,:);
reg_tbl_order = reg_tbl(strcmp(reg_tbl.condition,'nonviol') & strcmp(reg_tbl.sound,'C') ,:);

%% Run GLM: trial type
clear glm_output


% Identity model
u_t_mdl_identity = []; anova_identity = [];
reg_tbl_all.sound = categorical(reg_tbl_all.sound);
reg_tbl_all.sound = reordercats(reg_tbl_all.sound, {'G', 'D', 'F', 'C','A'});
u_t_mdl_identity = fitlm(reg_tbl_all, 'firing_rate ~ sound');
anova_identity = anova(u_t_mdl_identity);

% Positional model
u_t_mdl_position = []; anova_position = [];
reg_tbl_order.order_pos = categorical(reg_tbl_order.order_pos);
reg_tbl_order.order_pos = reordercats(reg_tbl_order.order_pos, {'position_2','position_3','position_4','position_5'});
u_t_mdl_position = fitlm(reg_tbl_order, 'firing_rate ~ order_pos');
anova_position = anova(u_t_mdl_position);

% Backward transition model
u_t_mdl_backward = []; anova_backward = [];
u_t_mdl_backward = fitlm(reg_tbl_order, 'firing_rate ~ backward_prob');
anova_backward = anova(u_t_mdl_backward);

% Forward transition model
u_t_mdl_forward = []; anova_forward = [];
u_t_mdl_forward = fitlm(reg_tbl_order, 'firing_rate ~ forward_prob');
anova_forward = anova(u_t_mdl_forward);


%% Output GLM

glm_sequence_out = table({u_t_mdl_identity}, {u_t_mdl_position}, {u_t_mdl_backward}, {u_t_mdl_forward},...
    'VariableNames',{'identity_mdl','position_mdl','backward_mdl','forward_mdl'});

anova_sequence_out = table({anova_identity}, {anova_position}, {anova_backward}, {anova_forward},...
    'VariableNames',{'anova_identity','anova_position','anova_backward','anova_forward'});
