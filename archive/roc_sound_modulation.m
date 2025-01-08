function [encoding_flag, encoding_roc, z_sdf] = roc_sound_modulation(sound_info_in)

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
u_t_mdl = [];

valid_trial_idx = [];
valid_trial_idx = find(reg_tbl.valid == 1 );

reg_tbl_trialtype = [];
reg_tbl_trialtype = reg_tbl(valid_trial_idx ,:);

sound_list = {'A','C','D','F','G'};

clear roc_x fr_* roc_out p_out

% For each averaged time point
for timepoint_i = 1:n_times

    % Input the timepoint specific firing times
    reg_tbl_trialtype.firing_rate = window_sdf(valid_trial_idx,timepoint_i);

    fr_baseline = [];
    fr_baseline = reg_tbl_trialtype.firing_rate(strcmp(reg_tbl_trialtype.condition,'Baseline'));

    for sound_i = 1:5

        fr_sound = []; fr_sound = reg_tbl_trialtype.firing_rate(strcmp(reg_tbl_trialtype.sound,sound_list{sound_i}));

        clear roc_x
        try
            [roc_x] = roc_curve(fr_baseline, fr_sound);
            roc_out(sound_i,timepoint_i) = roc_x.param.AROC;
            p_out(sound_i,timepoint_i) = permutationTest(fr_baseline, fr_sound, 100);
        catch
            roc_out(sound_i,timepoint_i) = NaN;
            p_out(sound_i,timepoint_i) = 0;
        end

    end

    for pos_i = 1:5

        fr_pos = []; fr_pos = reg_tbl_trialtype.firing_rate(strcmp(reg_tbl_trialtype.order_pos,['position_' int2str(pos_i)]));

        clear roc_x

        try
            [roc_x] = roc_curve(fr_baseline, fr_pos);
            roc_out(pos_i+5,timepoint_i) = roc_x.param.AROC;
            p_out(pos_i+5,timepoint_i) = permutationTest(fr_baseline, fr_pos, 100);
        catch
            roc_out(pos_i+5,timepoint_i) = NaN;
            p_out(pos_i+5,timepoint_i) = 0;
        end
    end
end


%% Determine periods of significance

signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 50 ms of selectivity following element

% Trial-type
clear start sig_len
for cond_i = 1:10
    [start{cond_i}, sig_len{cond_i}, ~] = ZeroOnesCount(roc_out(cond_i,analysis_win_idx)); % choice direction
end

encoding_flag = []; encoding_roc = [];

for cond_i = 1:10
    encoding_flag(1,cond_i) = any(sig_len{cond_i} >= signal_detect_wins);

    if encoding_flag(1,cond_i) == 1
        sig_range = [];
        sig_range(1) = analysis_win_idx(1)+start{cond_i}(find(sig_len{cond_i} >= signal_detect_wins,1,'first'));
        sig_range(2) = sig_range(1)+sig_len{cond_i}(find(sig_len{cond_i} >= signal_detect_wins,1,'first'));


        encoding_roc(1,cond_i) = nanmean(roc_out(cond_i,sig_range(1):sig_range(2)));
    else
        encoding_roc(1,cond_i) = NaN;
    end
end


end