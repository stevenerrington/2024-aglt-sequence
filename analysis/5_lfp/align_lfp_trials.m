function [lfp_trial] = align_lfp_trials(signal_in,alignTimes,timeWin)


lfp_trial = nan(length(alignTimes),range(timeWin)+1);

for ii = 1:length(alignTimes)
    try
        lfp_trial(ii,:) = signal_in(alignTimes(ii)+timeWin(1):alignTimes(ii)+timeWin(end));
    end
end
