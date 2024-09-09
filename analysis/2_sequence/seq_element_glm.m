clear glm_seq* encoding_seq_*

% WE DONT NEED A MOVING WINDOW - JUST DO A FIXED WINDOW AND CALL IT A DAY
neuron_i = 1498;
sound_info_in = sdf_soundAlign_data{neuron_i};

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
        reg_tbl.forward_prob(trial_i) = transitional_probability.forward_prob(curr_trial_ele,prev_trial_ele);
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

%% Run GLM: trial type
clear glm_output

% Initialise different model tables/data
reg_tbl_all = reg_tbl(strcmp(reg_tbl.condition,'nonviol') ,:);
reg_tbl_order = reg_tbl(strcmp(reg_tbl.condition,'nonviol') & strcmp(reg_tbl.sound,'C') ,:);

% Identity model
u_t_mdl_identity = [];
reg_tbl_all.sound = categorical(reg_tbl_all.sound);
reg_tbl_all.sound = reordercats(reg_tbl_all.sound, {'G', 'D', 'F', 'C','A'});
u_t_mdl_identity = fitlm(reg_tbl_all, 'firing_rate ~ exp_i + sound');

% Positional model
u_t_mdl_position = [];
reg_tbl_order.order_pos = categorical(reg_tbl_order.order_pos);
reg_tbl_order.order_pos = reordercats(reg_tbl_order.order_pos, {'position_2','position_3','position_4','position_5'});
u_t_mdl_position = fitlm(reg_tbl_order, 'firing_rate ~ exp_i + order_pos');

% Backward transition model
u_t_mdl_backward = [];
%! u_t_mdl_backward = fitlm(reg_tbl_order, 'firing_rate ~ exp_i + backward_prob');

% Forward transition model
u_t_mdl_forward = [];
%! u_t_mdl_forward = fitlm(reg_tbl_order, 'firing_rate ~ exp_i + forward_prob');




% GLM output -------------------------------
% - Ordinal Position
glm_output.ordinal.sig_times(1,timepoint_i) = u_t_mdl_position.Coefficients.pValue(2) < .01; % Pos 3
glm_output.ordinal.beta_weights(1,timepoint_i) = u_t_mdl_position.Coefficients.tStat(2); % Pos 3
glm_output.relative.var_exp(1,timepoint_i) = u_t_mdl_relative.Rsquared.Ordinary;
glm_output.relative.LogLikelihood(1,timepoint_i) = u_t_mdl_relative.LogLikelihood;

