function [sig_flag, dir_flag] = element_detect_modZ(normFR_in, spike_log)


norm_fr_soundA = normFR_in.norm_fr_soundA;
norm_fr_soundC = normFR_in.norm_fr_soundC;
norm_fr_soundD = normFR_in.norm_fr_soundD;
norm_fr_soundF = normFR_in.norm_fr_soundF;
norm_fr_soundG = normFR_in.norm_fr_soundG;

%% Analysis: determine significant modulation
% Define significance threshold and minimum significant time duration
sig_threshold = 5;
ops.min_sig_time = 50;

% Clear any previous significance and direction flags
clear sig_flag dir_flag

% Loop through each neuron in the spike log
for neuron_i = 1:size(spike_log,1)

    % Identify significant time points where normalized firing rate exceeds the threshold
    sig_idx_a = abs(norm_fr_soundA(neuron_i,200+[0:400])) >= sig_threshold;
    sig_idx_c = abs(norm_fr_soundC(neuron_i,200+[0:400])) >= sig_threshold;
    sig_idx_d = abs(norm_fr_soundD(neuron_i,200+[0:400])) >= sig_threshold;
    sig_idx_f = abs(norm_fr_soundF(neuron_i,200+[0:400])) >= sig_threshold;
    sig_idx_g = abs(norm_fr_soundG(neuron_i,200+[0:400])) >= sig_threshold;

    % Calculate the mean firing rate for each sound type
    fr_idx_a = nanmean(norm_fr_soundA(neuron_i,200+[0:400]));
    fr_idx_c = nanmean(norm_fr_soundC(neuron_i,200+[0:400]));
    fr_idx_d = nanmean(norm_fr_soundD(neuron_i,200+[0:400]));
    fr_idx_f = nanmean(norm_fr_soundF(neuron_i,200+[0:400]));
    fr_idx_g = nanmean(norm_fr_soundG(neuron_i,200+[0:400]));

    % Use ZeroOnesCount function to find the start and length of consecutive significant periods
    [start_a, len_a] = ZeroOnesCount(sig_idx_a);
    [start_c, len_c] = ZeroOnesCount(sig_idx_c);
    [start_d, len_d] = ZeroOnesCount(sig_idx_d);
    [start_f, len_f] = ZeroOnesCount(sig_idx_f);
    [start_g, len_g] = ZeroOnesCount(sig_idx_g);

    % Determine if there are any significant periods longer than the minimum time threshold
    idx_w_sigTimes_a = any(len_a > ops.min_sig_time);
    idx_w_sigTimes_c = any(len_c > ops.min_sig_time);
    idx_w_sigTimes_d = any(len_d > ops.min_sig_time);
    idx_w_sigTimes_f = any(len_f > ops.min_sig_time);
    idx_w_sigTimes_g = any(len_g > ops.min_sig_time);

    % Store the significance flags for each sound type
    sig_flag(neuron_i,:) = [idx_w_sigTimes_a, idx_w_sigTimes_c, idx_w_sigTimes_d, idx_w_sigTimes_f, idx_w_sigTimes_g];

    % Determine the direction of the firing rate (positive or negative)
    dir_flag(neuron_i,:) = [fr_idx_a > 0, fr_idx_c > 0, fr_idx_d > 0, fr_idx_f > 0, fr_idx_g > 0];
end


clear sig_idx* fr_* start_* len_* idx_w*