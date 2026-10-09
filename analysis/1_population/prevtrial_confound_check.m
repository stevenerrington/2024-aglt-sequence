%% Previous-trial carryover check: does trial N-1's sequence group predict trial N's?
%
% This is the decisive check for the "carryover/priming from the previous
% trial" explanation of the identical-stimulus confound
% (sequence_confound_check.m: frontal 89% / auditory 67% decoding of the
% CURRENT trial's sequence group from position-1 activity, where the
% stimulus is physically identical across groups).
%
% For a pure previous-trial echo in neural activity to explain THAT
% result, trial N-1's group must be statistically related to trial N's
% group -- otherwise an echo of trial N-1 in trial N's activity carries
% no information about trial N's group at all, and cannot account for
% decoding trial N's group specifically. So this checks the ONE thing
% that determines whether the carryover hypothesis is even viable, before
% spending time on a neural re-run:
%   - if trial N-1's group is independent of trial N's group: carryover
%     is ruled out as an explanation for the original finding.
%   - if trial N-1's group predicts trial N's group (e.g. a "no immediate
%     repeat" quota-generation rule creating negative autocorrelation):
%     carryover remains viable, and is worth testing directly in neural
%     activity next.
%
% No neural data needed -- only each session's event table, matched trial
% N to trial N-1 by literal trial_n (so the "previous trial" is whatever
% actually happened immediately before, including violation/error trials
% -- only pairs where BOTH trial N and trial N-1 were valid nonviolation
% trials with a defined sequence group are used).
%
% Assumes: dirs, spike_log already in the workspace.

seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};

unique_sessions = unique(spike_log.session);
n_sessions = length(unique_sessions);

prev_group_all = [];
curr_group_all = [];

for session_i = 1:n_sessions

    try
        event_table_in = load(fullfile(dirs.mat_data, [unique_sessions{session_i} '.mat']), 'event_table');
    catch
        continue
    end

    event_table = event_table_in.event_table;
    nonviol_mask = strcmp(event_table.cond_label,'nonviol') & ~isnan(event_table.rewardOnset_ms);

    trial_n_all_session = event_table.trial_n;
    cond_value_session  = event_table.cond_value;

    group_session = nan(size(cond_value_session));
    for g = 1:4
        group_session(ismember(cond_value_session, seq_cond_values{g})) = g;
    end
    group_session(~nonviol_mask) = nan;   % only nonviolation trials have a defined group

    % Build a lookup: for each literal trial_n value in this session, what
    % group (if any) was it?
    max_trial_n = max(trial_n_all_session);
    group_by_trialn = nan(max_trial_n,1);
    group_by_trialn(trial_n_all_session) = group_session;

    valid_curr = find(~isnan(group_session));
    for k = 1:length(valid_curr)
        row_i = valid_curr(k);
        this_trial_n = trial_n_all_session(row_i);
        if this_trial_n <= 1
            continue
        end
        prev_group = group_by_trialn(this_trial_n - 1);
        if ~isnan(prev_group)
            prev_group_all = [prev_group_all; prev_group]; %#ok<AGROW>
            curr_group_all = [curr_group_all; group_session(row_i)]; %#ok<AGROW>
        end
    end

end

fprintf('Consecutive nonviolation trial pairs (trial N-1 -> trial N), both with a defined group: %d\n', ...
    length(curr_group_all));

%% Contingency table: previous group (rows) x current group (columns)

contingency = zeros(4,4);
for a = 1:4
    for b = 1:4
        contingency(a,b) = sum(prev_group_all == a & curr_group_all == b);
    end
end

fprintf('\nPrevious-trial group (rows) x current-trial group (columns):\n');
disp(contingency)

[chi2_stat, df, p_value, expected] = chi2_independence(contingency);
fprintf('chi2(%d) = %.2f, p = %.4g\n', df, chi2_stat, p_value);

% Same-group repeat rate, for an intuitive read alongside the formal test
repeat_rate = mean(prev_group_all == curr_group_all);
fprintf('Observed same-group repeat rate: %.3f (chance if independent = 0.250)\n', repeat_rate);

if p_value < 0.05
    fprintf(['\n--> Trial N-1''s group DOES predict trial N''s group.\n' ...
             '    Carryover from the previous trial remains a viable explanation for the\n' ...
             '    identical-stimulus confound -- worth testing directly in neural activity next.\n']);
else
    fprintf(['\n--> Trial N-1''s group does NOT predict trial N''s group (independent).\n' ...
             '    A previous-trial echo in neural activity, however strong, could not by itself\n' ...
             '    explain decoding of the CURRENT trial''s group -- this rules out the carryover\n' ...
             '    explanation for the original finding.\n']);
end
